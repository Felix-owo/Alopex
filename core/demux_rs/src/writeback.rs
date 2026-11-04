//! Stage spool与最终FASTQ共用的有界回写、句柄所有权及失败传播。
use crate::{Arc, Context, DemuxMode, JoinHandle, PathBuf, Result};
use rayon::prelude::*;
use std::collections::VecDeque;
use std::fs::File;
use std::io::Write;
#[cfg(test)]
use std::io::BufWriter;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Mutex, Weak};
use std::time::Duration;

pub(super) struct WritebackState {
    files: Mutex<VecDeque<Weak<WritebackFile>>>,
    max_entries: usize,
    pub(super) pressure: Option<(PathBuf, u64)>,
    stopped: AtomicBool,
    failed: AtomicBool,
    failure: Mutex<Option<String>>,
    written: AtomicU64,
    settled: AtomicU64,
    cycles: AtomicU64,
    synced_files: AtomicU64,
    sync_nanos: AtomicU64,
    peak_dirty: AtomicU64,
}

struct WritebackFile {
    file: File,
    written: AtomicU64,
    synced: AtomicU64,
    state: Weak<WritebackState>,
}

impl Drop for WritebackFile {
    fn drop(&mut self) {
        if let Some(state) = self.state.upgrade() {
            state.settled.fetch_add(
                self.written
                    .load(Ordering::Acquire)
                    .saturating_sub(self.synced.load(Ordering::Acquire)),
                Ordering::Release,
            );
        }
    }
}

pub(super) struct OutputFile {
    target: Arc<WritebackFile>,
    state: Option<Arc<WritebackState>>,
}

impl OutputFile {
    pub(super) fn new(file: File, state: Option<Arc<WritebackState>>) -> Result<Self> {
        let target = Arc::new(WritebackFile {
            file,
            written: AtomicU64::new(0),
            synced: AtomicU64::new(0),
            state: state.as_ref().map(Arc::downgrade).unwrap_or_default(),
        });
        if let Some(state) = &state {
            let mut files = state
                .files
                .lock()
                .map_err(|_| anyhow::anyhow!("writeback registry poisoned"))?;
            if files.len() >= state.max_entries {
                files.retain(|entry| entry.strong_count() > 0);
            }
            files.push_back(Arc::downgrade(&target));
        }
        Ok(Self { target, state })
    }
}

pub(super) fn spool_output(file: File, guard: Option<&WritebackGuard>) -> Result<OutputFile> {
    OutputFile::new(file, guard.map(|guard| guard.state.clone()))
}

impl OutputFile {
    pub(super) fn sync_spool(&self) -> std::io::Result<()> {
        self.target.file.sync_all()?;
        release_output_cache(&self.target.file);
        Ok(())
    }
}

impl Write for OutputFile {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if let Some(state) = &self.state {
            if state.failed.load(Ordering::Acquire) {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    state.error_message(),
                ));
            }
        }
        let count = (&self.target.file).write(bytes)?;
        if let Some(state) = &self.state {
            self.target
                .written
                .fetch_add(count as u64, Ordering::Release);
            state.written.fetch_add(count as u64, Ordering::Release);
        }
        Ok(count)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        (&self.target.file).flush()
    }
}

fn memory_limit(text: &str) -> Option<u64> {
    text.trim().parse::<u64>().ok().filter(|value| *value > 0)
}

fn writeback_pressure() -> Option<(PathBuf, u64)> {
    let membership = std::fs::read_to_string("/proc/self/cgroup").ok()?;
    let relative = membership
        .lines()
        .find_map(|line| line.strip_prefix("0::"))?;
    let base = Path::new("/sys/fs/cgroup");
    let current = base.join(relative.trim_start_matches('/'));
    if current
        .components()
        .any(|part| matches!(part, std::path::Component::ParentDir))
    {
        return None;
    }
    let mut selected: Option<(PathBuf, u64)> = None;
    for directory in current
        .ancestors()
        .take_while(|path| path.starts_with(base))
    {
        for name in ["memory.max", "memory.high"] {
            let limit = std::fs::read_to_string(directory.join(name))
                .ok()
                .and_then(|text| memory_limit(&text));
            if let Some(limit) = limit {
                if selected.as_ref().map_or(true, |(_, old)| limit < *old) {
                    selected = Some((directory.join("memory.stat"), limit));
                }
            }
        }
    }
    selected.map(|(path, limit)| {
        (
            path,
            (limit / 4).clamp(64 * 1024 * 1024, 6 * 1024 * 1024 * 1024),
        )
    })
}

fn dirty_bytes(text: &str) -> Option<u64> {
    text.lines().find_map(|line| {
        let mut fields = line.split_whitespace();
        (fields.next()? == "file_dirty")
            .then(|| fields.next()?.parse().ok())
            .flatten()
    })
}

#[cfg(target_os = "linux")]
fn release_output_cache(file: &File) {
    use std::os::fd::AsRawFd;
    unsafe extern "C" {
        fn posix_fadvise(fd: i32, offset: i64, length: i64, advice: i32) -> i32;
    }
    unsafe {
        let _ = posix_fadvise(file.as_raw_fd(), 0, 0, 4);
    }
}

#[cfg(not(target_os = "linux"))]
fn release_output_cache(_file: &File) {}

impl WritebackState {
    fn new(max_streams: usize, pressure: Option<(PathBuf, u64)>) -> Self {
        Self {
            files: Mutex::new(VecDeque::new()),
            max_entries: max_streams.saturating_add(128).max(256),
            pressure,
            stopped: AtomicBool::new(false),
            failed: AtomicBool::new(false),
            failure: Mutex::new(None),
            written: AtomicU64::new(0),
            settled: AtomicU64::new(0),
            cycles: AtomicU64::new(0),
            synced_files: AtomicU64::new(0),
            sync_nanos: AtomicU64::new(0),
            peak_dirty: AtomicU64::new(0),
        }
    }

    fn error_message(&self) -> String {
        self.failure
            .lock()
            .ok()
            .and_then(|error| error.clone())
            .unwrap_or_else(|| "output writeback failed".to_string())
    }

    fn pressured(&self) -> bool {
        if let Some((path, threshold)) = &self.pressure {
            if let Some(dirty) = std::fs::read_to_string(path)
                .ok()
                .and_then(|text| dirty_bytes(&text))
            {
                self.peak_dirty.fetch_max(dirty, Ordering::Relaxed);
                if dirty >= *threshold {
                    return true;
                }
            }
        }
        self.written
            .load(Ordering::Acquire)
            .saturating_sub(self.settled.load(Ordering::Acquire))
            >= 512 * 1024 * 1024
    }

    fn sync_file(&self, entry: &Weak<WritebackFile>) -> Result<()> {
        let Some(target) = entry.upgrade() else {
            return Ok(());
        };
        let written = target.written.load(Ordering::Acquire);
        let previous = target.synced.load(Ordering::Acquire);
        if written == previous {
            return Ok(());
        }
        let at = std::time::Instant::now();
        target
            .file
            .sync_data()
            .context("Failed writing back Droplet FASTQ data")?;
        release_output_cache(&target.file);
        target.synced.store(written, Ordering::Release);
        self.settled
            .fetch_add(written.saturating_sub(previous), Ordering::Release);
        self.synced_files.fetch_add(1, Ordering::Relaxed);
        self.sync_nanos.fetch_add(
            at.elapsed().as_nanos().min(u64::MAX as u128) as u64,
            Ordering::Relaxed,
        );
        Ok(())
    }

    fn pending_bytes(entry: &Weak<WritebackFile>) -> u64 {
        entry
            .upgrade()
            .map(|target| {
                target
                    .written
                    .load(Ordering::Acquire)
                    .saturating_sub(target.synced.load(Ordering::Acquire))
            })
            .unwrap_or(0)
    }

    fn run(&self) -> Result<()> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .thread_name(|index| format!("droplet-wb-{index}"))
            .build()
            .context("Cannot initialize Droplet writeback workers")?;
        while !self.stopped.load(Ordering::Acquire) {
            if self.pressured() {
                let batch = {
                    let mut files = self
                        .files
                        .lock()
                        .map_err(|_| anyhow::anyhow!("writeback registry poisoned"))?;
                    let mut ranked = files
                        .drain(..)
                        .map(|entry| {
                            let pending = Self::pending_bytes(&entry);
                            (entry, pending)
                        })
                        .filter(|(entry, _)| entry.strong_count() > 0)
                        .collect::<Vec<_>>();
                    ranked.sort_unstable_by_key(|(_, pending)| std::cmp::Reverse(*pending));
                    let count = ranked
                        .partition_point(|(_, pending)| *pending > 0)
                        .min(1024);
                    let rest = ranked.split_off(count);
                    files.extend(rest.into_iter().map(|(entry, _)| entry));
                    ranked
                        .into_iter()
                        .map(|(entry, _)| entry)
                        .collect::<Vec<_>>()
                };
                pool.install(|| batch.par_iter().try_for_each(|entry| self.sync_file(entry)))?;
                let mut files = self
                    .files
                    .lock()
                    .map_err(|_| anyhow::anyhow!("writeback registry poisoned"))?;
                files.extend(batch.into_iter().filter(|entry| entry.strong_count() > 0));
                self.cycles.fetch_add(1, Ordering::Relaxed);
            }
            std::thread::park_timeout(Duration::from_millis(250));
        }
        Ok(())
    }
}

pub(super) struct WritebackGuard {
    pub(super) state: Arc<WritebackState>,
    thread: Option<JoinHandle<Result<()>>>,
}

impl WritebackGuard {
    pub(super) fn start(mode: DemuxMode, max_streams: usize) -> Result<Option<Self>> {
        if mode != DemuxMode::DnaOnlyDroplet || !cfg!(target_os = "linux") {
            return Ok(None);
        }
        Self::spawn(Arc::new(WritebackState::new(
            max_streams,
            writeback_pressure(),
        )))
        .map(Some)
    }

    fn spawn(state: Arc<WritebackState>) -> Result<Self> {
        let worker = state.clone();
        let thread = std::thread::Builder::new()
            .name("droplet-writeback".to_string())
            .spawn(move || {
                let result =
                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| worker.run()))
                        .unwrap_or_else(|_| {
                            Err(anyhow::anyhow!("Droplet writeback worker panicked"))
                        });
                if let Err(error) = &result {
                    if let Ok(mut stored) = worker.failure.lock() {
                        *stored = Some(format!("{error:#}"));
                    }
                    worker.failed.store(true, Ordering::Release);
                }
                result
            })
            .context("Cannot start Droplet writeback coordinator")?;
        Ok(Self {
            state,
            thread: Some(thread),
        })
    }

    pub(super) fn finish(mut self) -> Result<()> {
        self.state.stopped.store(true, Ordering::Release);
        if let Some(thread) = self.thread.take() {
            thread.thread().unpark();
            thread
                .join()
                .map_err(|_| anyhow::anyhow!("Droplet writeback coordinator panicked"))??;
        }
        eprintln!(
            "DROPLET_WRITEBACK_PROFILE {}",
            serde_json::json!({
                "cycles":self.state.cycles.load(Ordering::Relaxed),
                "files_synced":self.state.synced_files.load(Ordering::Relaxed),
                "sync_seconds_sum":self.state.sync_nanos.load(Ordering::Relaxed) as f64 / 1e9,
                "peak_observed_dirty_bytes":self.state.peak_dirty.load(Ordering::Relaxed),
                "dirty_trigger_bytes":self.state.pressure.as_ref().map(|(_, value)| value),
                "timing_note":"background intervals overlap scientific processing"
            })
        );
        Ok(())
    }
}

impl Drop for WritebackGuard {
    fn drop(&mut self) {
        self.state.stopped.store(true, Ordering::Release);
        if let Some(thread) = self.thread.take() {
            thread.thread().unpark();
            let _ = thread.join();
        }
    }
}

#[test]
fn writeback_preserves_bytes_and_releases_weak_handles() {
    let root = std::env::temp_dir().join(format!("demux-writeback-bytes-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let path = root.join("output.gz");
    let state = Arc::new(WritebackState::new(8, None));
    let mut output = BufWriter::with_capacity(
        8,
        OutputFile::new(File::create(&path).unwrap(), Some(state.clone())).unwrap(),
    );
    output.write_all(b"first gzip member").unwrap();
    output.write_all(b"second gzip member").unwrap();
    output.flush().unwrap();
    let entry = state.files.lock().unwrap().front().unwrap().clone();
    state.sync_file(&entry).unwrap();
    assert_eq!(
        std::fs::read(&path).unwrap(),
        b"first gzip membersecond gzip member"
    );
    assert_eq!(
        state.written.load(Ordering::Acquire),
        state.settled.load(Ordering::Acquire)
    );
    output.write_all(b"tail").unwrap();
    output.flush().unwrap();
    drop(output);
    assert!(entry.upgrade().is_none());
    assert_eq!(
        state.written.load(Ordering::Acquire),
        state.settled.load(Ordering::Acquire)
    );
    assert_eq!(
        std::fs::read(&path).unwrap(),
        b"first gzip membersecond gzip membertail"
    );
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn writeback_registry_does_not_retain_evicted_output_fds() {
    let root = std::env::temp_dir().join(format!("demux-writeback-weak-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let path = root.join("output.gz");
    let state = Arc::new(WritebackState::new(8, None));
    for _ in 0..600 {
        let output = OutputFile::new(File::create(&path).unwrap(), Some(state.clone())).unwrap();
        assert_eq!(Arc::strong_count(&output.target), 1);
        drop(output);
    }
    assert!(state.files.lock().unwrap().len() <= state.max_entries);
    assert!(state
        .files
        .lock()
        .unwrap()
        .iter()
        .all(|entry| entry.upgrade().is_none()));
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn writeback_pressure_and_failure_propagate_before_completion() {
    assert_eq!(memory_limit("max\n"), None);
    assert_eq!(memory_limit("17179869184\n"), Some(16 * 1024 * 1024 * 1024));
    assert_eq!(
        dirty_bytes("anon 20\nfile_dirty 34\nfile_writeback 9\n"),
        Some(34)
    );
    assert_eq!(dirty_bytes("anon 20\n"), None);
    let root = std::env::temp_dir().join(format!("demux-writeback-failure-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let pressure = root.join("memory.stat");
    std::fs::write(&pressure, b"file_dirty 100\n").unwrap();
    let low_dirty = Arc::new(WritebackState::new(8, Some((pressure.clone(), 1000))));
    assert!(!low_dirty.pressured());
    low_dirty
        .written
        .store(512 * 1024 * 1024, Ordering::Release);
    assert!(low_dirty.pressured());
    low_dirty
        .settled
        .store(512 * 1024 * 1024, Ordering::Release);
    assert!(!low_dirty.pressured());
    let state = Arc::new(WritebackState::new(8, Some((pressure, 1))));
    let mut sink = OutputFile::new(File::open("/dev/null").unwrap(), Some(state.clone())).unwrap();
    sink.target.written.store(1, Ordering::Release);
    state.written.store(1, Ordering::Release);
    let guard = WritebackGuard::spawn(state.clone()).unwrap();
    let deadline = std::time::Instant::now() + Duration::from_secs(5);
    while !state.failed.load(Ordering::Acquire) && std::time::Instant::now() < deadline {
        std::thread::sleep(Duration::from_millis(10));
    }
    let error = guard.finish().unwrap_err();
    assert!(format!("{error:#}").contains("Failed writing back Droplet FASTQ data"));
    assert!(sink.write_all(b"must not be accepted").is_err());
    std::fs::remove_dir_all(root).unwrap();
}

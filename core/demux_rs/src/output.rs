use crate::writeback::{OutputFile, WritebackGuard, WritebackState};
type FastqWriter = BufWriter<OutputFile>;
use crate::compression::compress_plain_buffer;

use super::*;
use std::collections::VecDeque;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::process::Command;

pub(super) fn compress_accumulation(acc: &mut AccumulationBuffer) -> Result<()> {
    let mut compressor =
        Compressor::new(CompressionLvl::new(1).expect("compression level 1 is supported"));
    for (key, bucket) in acc.map.iter_mut() {
        if bucket.count == 0 || matches!(key, BucketKey::Unmatched) {
            continue;
        }
        compress_buffer(
            &mut compressor,
            &bucket.r1_raw,
            &mut bucket.r1_gzip,
            key,
            "R1",
        )?;
        compress_buffer(
            &mut compressor,
            &bucket.r2_raw,
            &mut bucket.r2_gzip,
            key,
            "R2",
        )?;
    }
    Ok(())
}

fn compress_buffer(
    compressor: &mut Compressor,
    input: &[u8],
    output: &mut Vec<u8>,
    key: &BucketKey,
    mate: &str,
) -> Result<()> {
    compress_plain_buffer(
        compressor,
        input,
        output,
        &format!("{} {}", key.label(), mate),
    )
}

pub(super) fn spawn_writer(
    rx_write: Receiver<CompletedChunk>,
    tx_acc_recycle: Sender<AccumulationBuffer>,
    paths: OutputPaths,
    cell_rows: Vec<CellRow>,
    completion: Receiver<bool>,
) -> JoinHandle<Result<u64>> {
    std::thread::spawn(move || {
        let buckets = paths.dna_cell_ids.len() + paths.rna_cell_ids.len();
        let writeback = WritebackGuard::start(
            paths.mode,
            writer_pair_budget(paths.mode, buckets, detect_soft_open_file_limit())
                .saturating_mul(2),
        )?;
        if writer_shards(paths.mode) > 1
            && writer_pair_budget(paths.mode, buckets, detect_soft_open_file_limit()) >= buckets
        {
            writer_loop_parallel(
                rx_write,
                tx_acc_recycle,
                paths,
                cell_rows,
                writeback,
                completion,
            )
        } else {
            writer_loop(
                rx_write,
                tx_acc_recycle,
                paths,
                cell_rows,
                writeback,
                completion,
            )
        }
    })
}

struct WriterCache {
    writers: AHashMap<BucketKey, (FastqWriter, FastqWriter)>,
    writeback: Option<Arc<WritebackState>>,
    order: VecDeque<BucketKey>,
    max_pairs: usize,
    all_buckets_fit: bool,
}

const MAX_OPEN_BUCKETS: usize = 128;
const WRITER_FD_RESERVE: usize = 64;

fn detect_soft_open_file_limit() -> Option<usize> {
    #[cfg(unix)]
    {
        let output = Command::new("/bin/sh")
            .arg("-c")
            .arg("ulimit -n")
            .output()
            .ok()?;
        if !output.status.success() {
            return None;
        }
        let text = String::from_utf8(output.stdout).ok()?;
        text.trim().parse::<usize>().ok()
    }
    #[cfg(not(unix))]
    {
        None
    }
}

fn safe_writer_pair_limit(requested: usize, soft_limit: Option<usize>) -> usize {
    let requested = requested.max(8);
    let Some(limit) = soft_limit else {
        return requested;
    };
    let safe_pairs = limit.saturating_sub(WRITER_FD_RESERVE) / 2;
    requested.min(safe_pairs.max(8))
}

fn writer_pair_budget(mode: DemuxMode, buckets: usize, soft_limit: Option<usize>) -> usize {
    let requested = if mode == DemuxMode::DnaOnlyDroplet && soft_limit.is_some() {
        buckets.max(MAX_OPEN_BUCKETS)
    } else {
        MAX_OPEN_BUCKETS
    };
    let reader_budget = if mode == DemuxMode::DnaOnlyDroplet {
        64
    } else {
        0
    };
    let available = safe_writer_pair_limit(
        requested,
        soft_limit.map(|limit| limit.saturating_sub(reader_budget)),
    );
    if available < buckets {
        available.min(MAX_OPEN_BUCKETS)
    } else {
        available
    }
}

#[test]
fn droplet_writer_budget_retains_all_or_uses_bounded_fallback() {
    let mode = DemuxMode::DnaOnlyDroplet;
    assert_eq!(writer_pair_budget(mode, 6000, Some(12128)), 6000);
    assert_eq!(writer_pair_budget(mode, 6000, Some(12127)), 128);
    assert_eq!(writer_pair_budget(mode, 6000, Some(256)), 64);
    assert_eq!(writer_pair_budget(mode, 6000, None), 128);
    assert_eq!(writer_pair_budget(mode, 13203, Some(65536)), 13203);
    assert_eq!(writer_pair_budget(mode, 40000, Some(65536)), 128);
    assert_eq!(writer_pair_budget(mode, 9000, Some(100000)), 9000);
    assert_eq!(
        writer_pair_budget(DemuxMode::DnaRna, 6000, Some(100000)),
        128
    );
}

impl WriterCache {
    fn new(mode: DemuxMode, buckets: usize, writeback: Option<Arc<WritebackState>>) -> Self {
        let soft_limit = detect_soft_open_file_limit();
        let max_pairs = MAX_OPEN_BUCKETS;
        let effective_max_pairs = writer_pair_budget(mode, buckets, soft_limit);
        if effective_max_pairs < max_pairs {
            log::warn!(
                "Reducing max open barcode buckets from {} to {} because the process open-file limit is {:?}; each bucket uses two FASTQ file descriptors",
                max_pairs,
                effective_max_pairs,
                soft_limit
            );
        }
        log::info!(
            "Writer cache: {} pairs for {} buckets; soft open-file limit {:?}",
            effective_max_pairs,
            buckets,
            soft_limit
        );
        Self {
            writers: AHashMap::new(),
            writeback,
            order: VecDeque::new(),
            max_pairs: effective_max_pairs,
            all_buckets_fit: mode == DemuxMode::DnaOnlyDroplet && buckets <= effective_max_pairs,
        }
    }

    fn touch(&mut self, key: &BucketKey) {
        if self.all_buckets_fit {
            return;
        }
        if let Some(pos) = self.order.iter().position(|existing| existing == key) {
            self.order.remove(pos);
        }
        self.order.push_back(key.clone());
    }

    fn get_or_open(
        &mut self,
        key: &BucketKey,
        paths: &OutputPaths,
    ) -> Result<&mut (FastqWriter, FastqWriter)> {
        if !self.writers.contains_key(key) {
            while self.writers.len() >= self.max_pairs {
                let old_key = self
                    .order
                    .pop_front()
                    .context("Writer cache order is empty while cache is full")?;
                if let Some((mut r1, mut r2)) = self.writers.remove(&old_key) {
                    r1.flush().with_context(|| {
                        format!("Failed flushing evicted R1 writer for {}", old_key.label())
                    })?;
                    r2.flush().with_context(|| {
                        format!("Failed flushing evicted R2 writer for {}", old_key.label())
                    })?;
                }
            }
            self.writers.insert(
                key.clone(),
                create_fastq_writers(key, paths, self.writeback.clone())?,
            );
        }
        self.touch(key);
        self.writers
            .get_mut(key)
            .with_context(|| format!("Writer missing after initialization for {}", key.label()))
    }

    fn flush_all(mut self) -> Result<()> {
        for (_, (mut r1, mut r2)) in self.writers.drain() {
            r1.flush().context("Failed flushing R1 writer")?;
            r2.flush().context("Failed flushing R2 writer")?;
        }
        Ok(())
    }
}

fn writer_loop(
    rx_write: Receiver<CompletedChunk>,
    tx_acc_recycle: Sender<AccumulationBuffer>,
    paths: OutputPaths,
    cell_rows: Vec<CellRow>,
    writeback: Option<WritebackGuard>,
    completion: Receiver<bool>,
) -> Result<u64> {
    let mut writers = WriterCache::new(
        paths.mode,
        paths.dna_cell_ids.len() + paths.rna_cell_ids.len(),
        writeback.as_ref().map(|guard| guard.state.clone()),
    );
    let mut stats: AHashMap<BucketKey, u64> = AHashMap::new();
    let mut matched_reads = 0u64;
    let mut fate_counts = FateCounts::default();

    let mut pending = BTreeMap::new();
    let mut next_chunk_id = 0u64;
    let profile_enabled = paths.mode == DemuxMode::DnaOnlyDroplet;
    let mut write_seconds = 0f64;
    let mut output_chunks = 0u64;
    for (chunk_id, buffer) in rx_write {
        for mut buffer in
            insert_completed_chunk(&mut pending, &mut next_chunk_id, chunk_id, buffer)?
        {
            let write_start = std::time::Instant::now();
            fate_counts.add_assign(buffer.fates);
            for (key, bucket) in buffer.map.iter_mut() {
                if bucket.count == 0 {
                    continue;
                }

                *stats.entry(key.clone()).or_insert(0) += bucket.count;
                if key.is_matched() {
                    matched_reads += bucket.count;
                } else {
                    continue;
                }

                let writer_pair = writers.get_or_open(key, &paths)?;
                writer_pair
                    .0
                    .write_all(&bucket.r1_gzip)
                    .with_context(|| format!("Failed writing R1 output for {}", key.label()))?;
                writer_pair
                    .1
                    .write_all(&bucket.r2_gzip)
                    .with_context(|| format!("Failed writing R2 output for {}", key.label()))?;
            }

            write_seconds += write_start.elapsed().as_secs_f64();
            output_chunks += 1;
            buffer.clear();
            tx_acc_recycle
                .send(buffer)
                .context("Accumulation buffer recycle channel is closed")?;
        }
    }
    if let Some(first_pending) = pending.keys().next() {
        anyhow::bail!(
            "Writer did not receive chunk {next_chunk_id}; first pending chunk is {first_pending}"
        );
    }

    let flush_start = std::time::Instant::now();
    writers.flush_all()?;
    if let Some(writeback) = writeback {
        writeback.finish()?;
    }
    if profile_enabled {
        eprintln!(
            "DROPLET_WRITER_PROFILE {}",
            serde_json::json!({
                "chunks":output_chunks, "write_and_accounting_seconds":write_seconds,
                "final_flush_seconds":flush_start.elapsed().as_secs_f64()
            })
        );
    }

    let fate_matched = fate_counts.matched_read_pairs();
    if fate_matched != matched_reads {
        anyhow::bail!(
            "Internal demux accounting mismatch: bucket matched_reads={} but fate matched_reads={}",
            matched_reads,
            fate_matched
        );
    }
    anyhow::ensure!(
        completion
            .recv()
            .context("Missing demux input completion")?,
        "Demux input processing failed"
    );
    write_report(&paths, &cell_rows, &stats, fate_counts)?;
    Ok(matched_reads)
}

fn writer_loop_parallel(
    rx: Receiver<CompletedChunk>,
    recycle: Sender<AccumulationBuffer>,
    paths: OutputPaths,
    rows: Vec<CellRow>,
    writeback: Option<WritebackGuard>,
    completion: Receiver<bool>,
) -> Result<u64> {
    let partitions = writer_shards(paths.mode);
    let mut senders = Vec::new();
    let mut workers = Vec::new();
    for partition in 0..partitions {
        let (tx, rx) = bounded::<(u64, Vec<(BucketKey, FastqBucket)>)>(2);
        senders.push(tx);
        let paths = paths.clone();
        let buckets = rows.len();
        let writeback_state = writeback.as_ref().map(|guard| guard.state.clone());
        workers.push(std::thread::spawn(move|| -> Result<()> {
                let mut cache=WriterCache::new(paths.mode,buckets,writeback_state);let mut next=0u64;let mut elapsed=0.0;
                for (id,groups) in rx {
                    anyhow::ensure!(id==next,"writer partition sequence mismatch");next+=1;let at=std::time::Instant::now();
                    for (key,bucket) in groups {
                        let pair=cache.get_or_open(&key,&paths)?;
                        pair.0.write_all(&bucket.r1_gzip).with_context(||format!("R1 partition write failed: {}",key.label()))?;
                        pair.1.write_all(&bucket.r2_gzip).with_context(||format!("R2 partition write failed: {}",key.label()))?;
                    }
                    elapsed+=at.elapsed().as_secs_f64();
                }
                let at=std::time::Instant::now();cache.flush_all()?;
                eprintln!("DROPLET_WRITER_PARTITION_PROFILE {}",serde_json::json!({"partition":partition,"blocks":next,"write_seconds":elapsed,"flush_seconds":at.elapsed().as_secs_f64()}));Ok(())
            }));
    }
    let mut stats: AHashMap<BucketKey, u64> = AHashMap::new();
    let mut fates = FateCounts::default();
    let mut matched = 0u64;
    let mut pending = BTreeMap::new();
    let mut next = 0u64;
    let mut dispatch = 0u64;
    let mut queue_seconds = 0.0;
    let processing = (|| -> Result<()> {
        for (id, buffer) in rx {
            for mut buffer in insert_completed_chunk(&mut pending, &mut next, id, buffer)? {
                fates.add_assign(buffer.fates);
                let mut groups: Vec<Vec<(BucketKey, FastqBucket)>> =
                    (0..partitions).map(|_| Vec::new()).collect();
                for (key, bucket) in buffer.map.drain() {
                    if bucket.count == 0 {
                        continue;
                    }
                    *stats.entry(key.clone()).or_default() += bucket.count;
                    if !key.is_matched() {
                        continue;
                    }
                    matched += bucket.count;
                    let mut hash = 0xcbf29ce484222325u64;
                    for byte in key.label().bytes() {
                        hash = (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3);
                    }
                    groups[(hash as usize) % partitions].push((key, bucket));
                }
                let at = std::time::Instant::now();
                for (tx, groups) in senders.iter().zip(groups) {
                    tx.send((dispatch, groups))?;
                }
                queue_seconds += at.elapsed().as_secs_f64();
                dispatch += 1;
                buffer.clear();
                recycle.send(buffer)?;
            }
        }
        anyhow::ensure!(pending.is_empty(), "missing parallel writer block");
        Ok(())
    })();
    drop(senders);
    let mut failure = None;
    for worker in workers {
        match worker.join() {
            Ok(Ok(())) => {}
            Ok(Err(error)) => {
                if failure.is_none() {
                    failure = Some(error);
                }
            }
            Err(_) => {
                if failure.is_none() {
                    failure = Some(anyhow::anyhow!("writer partition panic"));
                }
            }
        }
    }
    if let Some(error) = failure {
        return Err(error);
    }
    processing?;
    if let Some(writeback) = writeback {
        writeback.finish()?;
    }
    anyhow::ensure!(
        fates.matched_read_pairs() == matched,
        "parallel writer pair accounting mismatch"
    );
    eprintln!(
        "DROPLET_WRITER_PROFILE {}",
        serde_json::json!({"partitions":partitions,"chunks":dispatch,"dispatch_wait_seconds":queue_seconds})
    );
    anyhow::ensure!(
        completion
            .recv()
            .context("Missing demux input completion")?,
        "Demux input processing failed"
    );
    write_report(&paths, &rows, &stats, fates)?;
    Ok(matched)
}

fn droplet_writer_buffer_capacity(outputs: usize, dirty_threshold: Option<u64>) -> usize {
    let budget = dirty_threshold
        .map(|threshold| (threshold / 2).clamp(512 * 1024 * 1024, 2 * 1024 * 1024 * 1024))
        .unwrap_or(512 * 1024 * 1024);
    (budget as usize / outputs.max(1)).clamp(8 * 1024, 256 * 1024)
}

#[test]
fn droplet_buffer_budget_is_bounded_and_scales_for_large_cell_sets() {
    assert_eq!(
        droplet_writer_buffer_capacity(26406, None),
        512 * 1024 * 1024 / 26406
    );
    let large = droplet_writer_buffer_capacity(26406, Some(4 * 1024 * 1024 * 1024));
    assert_eq!(large, 2 * 1024 * 1024 * 1024 / 26406);
    assert!(large * 26406 <= 2 * 1024 * 1024 * 1024);
    assert_eq!(
        droplet_writer_buffer_capacity(2, Some(4 * 1024 * 1024 * 1024)),
        256 * 1024
    );
    assert_eq!(
        droplet_writer_buffer_capacity(1_000_000, Some(64 * 1024 * 1024)),
        8 * 1024
    );
}

fn create_fastq_writers(
    key: &BucketKey,
    paths: &OutputPaths,
    writeback: Option<Arc<WritebackState>>,
) -> Result<(FastqWriter, FastqWriter)> {
    let dir = key.output_dir(paths)?;
    std::fs::create_dir_all(&dir)
        .with_context(|| format!("Cannot create output directory: {}", dir.display()))?;
    let stem = key.file_stem(paths)?;
    let r1_path = dir.join(format!("{}_R1.fastq.gz", stem));
    let r2_path = dir.join(format!("{}_R2.fastq.gz", stem));
    let r1 = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&r1_path)
        .with_context(|| format!("Cannot open R1 output file: {}", r1_path.display()))?;
    let r2 = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&r2_path)
        .with_context(|| format!("Cannot open R2 output file: {}", r2_path.display()))?;
    let dirty_threshold = writeback
        .as_ref()
        .and_then(|state| state.pressure.as_ref())
        .map(|(_, limit)| *limit);
    let r1 = OutputFile::new(r1, writeback.clone())?;
    let r2 = OutputFile::new(r2, writeback)?;
    if paths.mode == DemuxMode::DnaOnlyDroplet {
        let outputs = (paths.dna_cell_ids.len() + paths.rna_cell_ids.len()).max(1) * 2;
        let capacity = droplet_writer_buffer_capacity(outputs, dirty_threshold);
        Ok((
            BufWriter::with_capacity(capacity, r1),
            BufWriter::with_capacity(capacity, r2),
        ))
    } else {
        Ok((BufWriter::new(r1), BufWriter::new(r2)))
    }
}

fn write_report(
    paths: &OutputPaths,
    cell_rows: &[CellRow],
    stats: &AHashMap<BucketKey, u64>,
    fate_counts: FateCounts,
) -> Result<()> {
    let mut rows = cell_rows.to_vec();
    rows.sort_by(|a, b| {
        match (a.cell_order.parse::<i64>(), b.cell_order.parse::<i64>()) {
            (Ok(left), Ok(right)) => left.cmp(&right),
            (Ok(_), Err(_)) => std::cmp::Ordering::Less,
            (Err(_), Ok(_)) => std::cmp::Ordering::Greater,
            (Err(_), Err(_)) => a.cell_order.cmp(&b.cell_order),
        }
        .then_with(|| a.plate_id.cmp(&b.plate_id))
    });

    let mut sample_stats = Vec::with_capacity(rows.len());
    for row in rows {
        let dna_count = *stats
            .get(&BucketKey::Dna(row.dna_barcode.clone()))
            .unwrap_or(&0);
        let rna_count = *stats
            .get(&BucketKey::Rna(row.rna_barcode.clone()))
            .unwrap_or(&0);
        let dna_status = output_status(paths, Modality::Dna, &row.dna_barcode, dna_count)?;
        let rna_status = if paths.mode == DemuxMode::DnaRna {
            output_status(paths, Modality::Rna, &row.rna_barcode, rna_count)?
        } else {
            "Not_Used".to_string()
        };
        let cell_sample_id = paths
            .dna_cell_ids
            .get(&row.dna_barcode)
            .with_context(|| format!("No canonical Sample ID for DNA barcode {}", row.dna_barcode))?
            .clone();
        sample_stats.push(SampleStats {
            sample_name: paths.sample_name.clone(),
            cell_sample_id,
            plate_id: row.plate_id.clone(),
            cell_order: row.cell_order.clone(),
            dna_barcode: row.dna_barcode.clone(),
            rna_barcode: row.rna_barcode.clone(),
            dna_read_count: dna_count,
            rna_read_count: rna_count,
            dna_status,
            rna_status,
        });
    }

    let usable_from_buckets: u64 = stats.values().sum();
    let usable_read_pairs = fate_counts.usable_read_pairs();
    if usable_from_buckets != usable_read_pairs {
        anyhow::bail!(
            "Internal demux accounting mismatch: bucket usable_reads={} but fate usable_reads={}",
            usable_from_buckets,
            usable_read_pairs
        );
    }
    let input_read_pairs = fate_counts.input_read_pairs();
    let matched_reads = fate_counts.matched_read_pairs();
    let report = DemuxReport {
        schema_version: 3,
        build_version: env!("CARGO_PKG_VERSION").to_string(),
        source_revision: BUILD_SOURCE_REVISION.to_string(),
        mode: match paths.mode {
            DemuxMode::DnaOnly => "dna-only".to_string(),
            DemuxMode::DnaOnlyTaps => "dna-only-taps".to_string(),
            DemuxMode::DnaOnlyDroplet => "dna-only-droplet".to_string(),
            DemuxMode::DnaRna => "dna-rna".to_string(),
        },
        retention_policy: "matched_read_pairs".to_string(),
        retention_threshold: paths.min_matched_read_pairs,
        input_fastq_pairs: paths.input_fastq_pairs,
        input_read_pairs,
        usable_read_pairs,
        matched_reads,
        read_fates: fate_counts,
        samples: sample_stats,
    };
    ensure_parent_dir(&paths.json_report, "--json-report")?;
    let file = File::create(&paths.json_report)
        .with_context(|| format!("Cannot create JSON report: {}", paths.json_report.display()))?;
    serde_json::to_writer_pretty(file, &report).context("Failed writing JSON report")?;
    Ok(())
}

#[test]
fn report_cell_order_has_transitive_numeric_then_text_order() {
    let root = std::env::temp_dir().join(format!(
        "alopex_report_order.{}.{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&root).unwrap();
    let expected = [
        "-9223372036854775808",
        "-10",
        "-2",
        "01",
        "1",
        "+2",
        "2",
        "10",
        "9223372036854775807",
        "1z",
        "A1",
        "B2",
    ];
    let rows: Vec<CellRow> = expected
        .iter()
        .enumerate()
        .map(|(i, order)| CellRow {
            dna_barcode: format!("barcode{i}"),
            rna_barcode: String::new(),
            plate_id: format!("P{i:02}"),
            cell_order: order.to_string(),
        })
        .collect();
    let paths = OutputPaths {
        dna_out_dir: root.join("dna"),
        rna_out_dir: root.join("rna"),
        json_report: root.join("report.json"),
        sample_name: "Sort".to_string(),
        dna_cell_ids: rows
            .iter()
            .map(|row| (row.dna_barcode.clone(), format!("Sort_{}", row.plate_id)))
            .collect(),
        rna_cell_ids: AHashMap::new(),
        min_matched_read_pairs: 1,
        mode: DemuxMode::DnaOnly,
        input_fastq_pairs: 1,
    };
    for reverse in [false, true] {
        for shift in 0..rows.len() {
            let mut input = rows.clone();
            input.rotate_left(shift);
            if reverse {
                input.reverse();
            }
            write_report(&paths, &input, &AHashMap::new(), FateCounts::default()).unwrap();
            let report: serde_json::Value =
                serde_json::from_slice(&std::fs::read(&paths.json_report).unwrap()).unwrap();
            let observed: Vec<&str> = report["samples"]
                .as_array()
                .unwrap()
                .iter()
                .map(|sample| sample["cell_order"].as_str().unwrap())
                .collect();
            assert_eq!(observed, expected, "reverse={reverse}, shift={shift}");
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}

fn output_status(
    paths: &OutputPaths,
    modality: Modality,
    barcode: &str,
    count: u64,
) -> Result<String> {
    if count == 0 {
        return Ok("Zero_Output".to_string());
    }
    let base_dir = match modality {
        Modality::Dna => &paths.dna_out_dir,
        Modality::Rna => &paths.rna_out_dir,
    };
    let key = match modality {
        Modality::Dna => BucketKey::Dna(barcode.to_string()),
        Modality::Rna => BucketKey::Rna(barcode.to_string()),
    };
    let cell_sample_id = key.cell_sample_id(paths)?;
    let dir = base_dir.join(cell_sample_id);
    let r1_path = dir.join(format!("{}_R1.fastq.gz", cell_sample_id));
    if !r1_path.exists() {
        return Ok("Zero_Output".to_string());
    }
    let below_threshold = count < paths.min_matched_read_pairs;
    if below_threshold {
        std::fs::remove_dir_all(&dir).with_context(|| {
            format!("Cannot remove low-read output directory: {}", dir.display())
        })?;
        return Ok("Low_Reads_Removed".to_string());
    }
    Ok("Pass".to_string())
}

pub(super) fn ensure_parent_dir(path: &Path, label: &str) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).with_context(|| {
                format!(
                    "Cannot create parent directory for {}: {}",
                    label,
                    parent.display()
                )
            })?;
        }
    }
    Ok(())
}

fn validate_cell_output_names(root: &Path, names: &[String]) -> Result<()> {
    let mut spellings: AHashMap<String, &String> = AHashMap::new();
    let mut candidates = AHashSet::new();
    for name in names {
        if let Some(previous) = spellings.insert(name.to_ascii_lowercase(), name) {
            candidates.insert(previous);
            candidates.insert(name);
        }
    }
    if candidates.is_empty() {
        return Ok(());
    }
    let paths: Vec<_> = candidates.into_iter().map(|name| root.join(name)).collect();
    let directories: Vec<_> = paths.iter().map(PathBuf::as_path).collect();
    crate::cli::validate_prospective_paths(&directories, &[], Some(root))
        .context("Cannot validate distinct cell output names")
}

pub(super) fn prepare_output_roots(config: &Config, spool: Option<&Path>, rows: &[CellRow]) -> Result<()> {
    let mut inputs: Vec<_> = config.r1.iter().map(|p| ("--r1", p.as_path()))
        .chain(config.r2.iter().map(|p| ("--r2", p.as_path())))
        .chain(std::iter::once(("--barcode-map", config.barcode_map.as_path())))
        .collect();
    if let Some(spool) = spool { inputs.push(("--droplet-spool", spool)); }
    let mut directories = vec![("--dna-out-dir", config.dna_out_dir.as_path())];
    if config.mode == DemuxMode::DnaRna {
        directories.push(("--rna-out-dir", config.rna_out_dir.as_path()));
    }
    crate::cli::validate_io_paths(&inputs, &directories,
        &[("--json-report", config.json_report.as_path())])?;
    let names: Vec<_> = rows.iter().map(|row| format!("{}_{}", config.sample_name,
        if config.mode == DemuxMode::DnaOnlyDroplet { &row.dna_barcode } else { &row.plate_id })).collect();
    for (_, root) in &directories {
        validate_cell_output_names(root, &names)?;
    }
    ensure_parent_dir(&config.json_report, "--json-report")?;
    for (label, path) in directories {
        let root = materialize_output_root(path, label)?;
        if root.exists() {
            let metadata = std::fs::symlink_metadata(&root)
                .with_context(|| format!("Cannot inspect {label}: {}", root.display()))?;
            if metadata.file_type().is_symlink() || !metadata.is_dir() {
                anyhow::bail!("{label} must be a real directory: {}", root.display());
            }
            std::fs::remove_dir_all(&root).with_context(|| {
                format!(
                    "Cannot reset exclusive run output root {label}: {}",
                    root.display()
                )
            })?;
        }
        std::fs::create_dir_all(&root)
            .with_context(|| format!("Cannot create {label}: {}", root.display()))?;
    }
    Ok(())
}

fn materialize_output_root(path: &Path, label: &str) -> Result<PathBuf> {
    if path.exists() {
        let metadata = std::fs::symlink_metadata(path)
            .with_context(|| format!("Cannot inspect {label}: {}", path.display()))?;
        if metadata.file_type().is_symlink() || !metadata.is_dir() {
            anyhow::bail!("{label} must be a real directory: {}", path.display());
        }
    } else {
        std::fs::create_dir_all(path)
            .with_context(|| format!("Cannot create {label}: {}", path.display()))?;
    }
    path.canonicalize()
        .with_context(|| format!("Cannot resolve {label}: {}", path.display()))
}

fn writer_shards(mode: DemuxMode) -> usize {
    if mode == DemuxMode::DnaOnlyDroplet {
        4
    } else {
        1
    }
}

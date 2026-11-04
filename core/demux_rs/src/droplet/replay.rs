//! 按原block/record顺序重放最终名单并提交配对输出；不参与calling。
use super::spool::*;
use crate::compression::compress_plain_buffer;
use crate::*;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::time::Instant;

fn compress_parallel(acc: &mut AccumulationBuffer, threads: usize) -> Result<()> {
    let mut buckets: Vec<_> = std::mem::take(&mut acc.map).into_iter().collect();
    let group = (buckets.len() + threads - 1) / threads;
    buckets.par_chunks_mut(group.max(1)).try_for_each(|group| {
        let mut compressor = Compressor::new(CompressionLvl::new(1).unwrap());
        for (key, bucket) in group {
            if bucket.count == 0 || !key.is_matched() {
                continue;
            }
            compress_plain_buffer(
                &mut compressor,
                &bucket.r1_raw,
                &mut bucket.r1_gzip,
                "R1 replay",
            )?;
            compress_plain_buffer(
                &mut compressor,
                &bucket.r2_raw,
                &mut bucket.r2_gzip,
                "R2 replay",
            )?;
            bucket.r1_raw = Vec::new();
            bucket.r2_raw = Vec::new();
        }
        Ok::<(), anyhow::Error>(())
    })?;
    acc.map = buckets.into_iter().collect();
    Ok(())
}
#[derive(Default, Serialize)]
struct ReplayTimings {
    decode_match_seconds: f64,
    assemble_seconds: f64,
    compression_seconds: f64,
    queue_wait_seconds: f64,
    blocks: u64,
}
fn consumed_spool_range(previous: u64, end: u64, force: bool) -> Option<(u64, u64)> {
    let through = end & !4095;
    let length = through.checked_sub(previous)?;
    if length == 0 || (!force && length < 8 * 1024 * 1024) {
        return None;
    }
    i64::try_from(previous).ok()?;
    i64::try_from(length).ok()?;
    Some((previous, length))
}

#[cfg(target_os = "linux")]
fn release_consumed_spool(file: &File, offset: u64, bytes: u64) -> bool {
    use std::os::fd::AsRawFd;
    unsafe extern "C" {
        fn posix_fadvise(fd: i32, offset: i64, length: i64, advice: i32) -> i32;
    }
    unsafe { posix_fadvise(file.as_raw_fd(), offset as i64, bytes as i64, 4) == 0 }
}

#[cfg(not(target_os = "linux"))]
fn release_consumed_spool(_file: &File, _offset: u64, _bytes: u64) -> bool {
    true
}

#[test]
fn consumed_spool_cache_ranges_wait_for_complete_batches_and_full_pages() {
    assert_eq!(consumed_spool_range(0, 4095, true), None);
    assert_eq!(consumed_spool_range(0, 8193, false), None);
    assert_eq!(consumed_spool_range(0, 8193, true), Some((0, 8192)));
    assert_eq!(consumed_spool_range(8192, 4096, true), None);
    assert_eq!(
        consumed_spool_range(0, 8 * 1024 * 1024 + 12, false),
        Some((0, 8 * 1024 * 1024))
    );
    assert_eq!(consumed_spool_range(u64::MAX - 4095, u64::MAX, true), None);
}

fn prepare_replay_block(
    block: &Block,
    readers: &[std::sync::Mutex<File>],
    indices: &Arc<Vec<BarcodeIndex>>,
    barcode_ids: &AHashMap<&str, usize>,
    barcode_list: &[String],
    threads: usize,
) -> Result<(u64, AccumulationBuffer, ReplayTimings)> {
    let mut timing = ReplayTimings::default();
    let start = Instant::now();
    let selected: Result<Vec<_>> = block
        .segments
        .par_iter()
        .enumerate()
        .map(|(shard, segment)| {
            let raw = {
                if segment.bytes == 0 {
                    anyhow::ensure!(segment.raw_bytes == 0, "empty compressed segment");
                    Vec::new()
                } else {
                    anyhow::ensure!(
                        segment.raw_bytes <= 1024 * 1024 * 1024
                            && segment.bytes <= 1024 * 1024 * 1024,
                        "oversized segment"
                    );
                    let mut gz = vec![0u8; segment.bytes as usize];
                    {
                        let mut file = readers[shard]
                            .lock()
                            .map_err(|_| anyhow::anyhow!("poisoned spool reader"))?;
                        file.seek(SeekFrom::Start(segment.offset))?;
                        file.read_exact(&mut gz)?;
                    }
                    let mut raw = vec![0u8; segment.raw_bytes as usize];
                    let n = libdeflate::Decompressor::new().gzip_decompress(&gz, &mut raw)?;
                    anyhow::ensure!(n == raw.len(), "decompressed size mismatch");
                    raw
                }
            };
            let mut selected = Vec::new();
            let mut fates = FateCounts::default();
            let mut seen = 0u64;
            walk(&raw, |record| {
                anyhow::ensure!(
                    u64::from(record.ordinal) < block.pairs,
                    "ordinal outside block"
                );
                seen += 1;
                match indices[0].match_simple(record.barcode) {
                    BarcodeMatch::Unique { id } if !record.short_r2 => {
                        fates.dna_assigned += 1;
                        let first = record.r1.as_ptr() as usize - raw.as_ptr() as usize;
                        let second = record.r2.as_ptr() as usize - raw.as_ptr() as usize;
                        selected.push((
                            record.ordinal,
                            barcode_ids[id],
                            shard,
                            first,
                            record.r1.len(),
                            second,
                            record.r2.len(),
                        ));
                    }
                    BarcodeMatch::Unique { .. } => fates.short += 1,
                    BarcodeMatch::Ambiguous => fates.ambiguous += 1,
                    BarcodeMatch::None => fates.unmatched += 1,
                }
                Ok(())
            })?;
            Ok((selected, fates, seen, raw))
        })
        .collect();
    let selected = selected?;
    timing.decode_match_seconds = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let mut acc = AccumulationBuffer::new();
    acc.fates.ambiguous = block.base_ambiguous;
    acc.fates.unmatched = block.base_unmatched;
    acc.fates.short = block.base_short;
    let mut records = Vec::new();
    let mut seen = block.base_ambiguous + block.base_unmatched + block.base_short;
    let mut raw_shards = Vec::new();
    for (mut selected, fates, n, raw) in selected {
        acc.fates.add_assign(fates);
        seen += n;
        records.append(&mut selected);
        raw_shards.push(raw);
    }
    anyhow::ensure!(
        seen == block.pairs && acc.fates.input_read_pairs() == block.pairs,
        "replay pair conservation"
    );
    records.sort_unstable_by_key(|r| r.0);
    anyhow::ensure!(
        records.windows(2).all(|w| w[0].0 < w[1].0),
        "duplicate selected ordinal"
    );
    let mut groups: Vec<Vec<_>> = (0..barcode_list.len()).map(|_| Vec::new()).collect();
    for record in records {
        groups[record.1].push(record);
    }
    let buckets: Vec<_> = groups
        .into_par_iter()
        .enumerate()
        .filter(|(_, group)| !group.is_empty())
        .map(|(bc, group)| {
            let mut bucket = FastqBucket::new();
            bucket.count = group.len() as u64;
            bucket.r1_raw.reserve(group.iter().map(|r| r.4).sum());
            bucket.r2_raw.reserve(group.iter().map(|r| r.6).sum());
            for (_, _, shard, first, n1, second, n2) in group {
                bucket
                    .r1_raw
                    .extend_from_slice(&raw_shards[shard][first..first + n1]);
                bucket
                    .r2_raw
                    .extend_from_slice(&raw_shards[shard][second..second + n2]);
            }
            (BucketKey::Dna(barcode_list[bc].clone()), bucket)
        })
        .collect();
    for (key, bucket) in buckets {
        acc.map.insert(key, bucket);
    }
    acc.bucket_mut(BucketKey::Unmatched).count = acc.fates.ambiguous + acc.fates.unmatched;
    timing.assemble_seconds = start.elapsed().as_secs_f64();
    let start = Instant::now();
    compress_parallel(&mut acc, threads)?;
    timing.compression_seconds = start.elapsed().as_secs_f64();
    timing.blocks = 1;
    Ok((block.id, acc, timing))
}
pub(crate) fn run(
    config: &Config,
    rows: Vec<CellRow>,
    indices: Arc<Vec<BarcodeIndex>>,
    paths: OutputPaths,
    spool: &Path,
) -> Result<()> {
    let (manifest, blocks) = blocks(spool)?;
    anyhow::ensure!(
        config.mode == DemuxMode::DnaOnlyDroplet
            && manifest.input_fastq_pairs == paths.input_fastq_pairs,
        "replay input mismatch"
    );
    let readers = (0..manifest.shards)
        .map(|shard| File::open(spool.join(format!("{shard:04}.pack"))).map(std::sync::Mutex::new))
        .collect::<std::io::Result<Vec<_>>>()?;
    let max_raw = blocks
        .iter()
        .map(|b| b.segments.iter().map(|s| s.raw_bytes).sum::<u64>())
        .max()
        .unwrap_or(1)
        .max(1);
    let window =
        (2usize * 1024 * 1024 * 1024 / (max_raw as usize).saturating_mul(4).max(1)).clamp(1, 4);
    let (tx, rx) = bounded(2);
    let (tx_recycle, rx_recycle) = bounded(2);
    let barcode_list: Vec<_> = rows.iter().map(|r| r.dna_barcode.clone()).collect();
    let barcode_ids: AHashMap<_, _> = barcode_list
        .iter()
        .enumerate()
        .map(|(i, bc)| (bc.as_str(), i))
        .collect();
    let drain = std::thread::spawn(move || {
        for buffer in rx_recycle {
            drop(buffer);
        }
    });
    let final_report = paths.json_report.clone();
    let private_report = final_report.with_extension("droplet_stage.partial");
    let paths = OutputPaths {
        json_report: private_report.clone(),
        ..paths
    };
    let (tx_completion, rx_completion) = bounded(1);
    let writer = spawn_writer(rx, tx_recycle, paths, rows, rx_completion);
    let start = Instant::now();
    let mut total = ReplayTimings::default();
    let mut released = vec![0u64; readers.len()];
    let mut cache_requested = 0u64;
    let mut cache_calls = 0u64;
    let mut cache_failed = 0u64;
    let processing = (|| -> Result<()> {
        for batch in blocks.chunks(window) {
            let prepared: Result<Vec<_>> = batch
                .par_iter()
                .map(|block| {
                    prepare_replay_block(
                        block,
                        &readers,
                        &indices,
                        &barcode_ids,
                        &barcode_list,
                        config.threads,
                    )
                })
                .collect();
            let prepared = prepared?;
            if cfg!(target_os = "linux") {
                let last = batch.last().unwrap();
                for (shard, segment) in last.segments.iter().enumerate() {
                    let Some(end) = segment.offset.checked_add(segment.bytes) else {
                        continue;
                    };
                    if let Some((offset, length)) =
                        consumed_spool_range(released[shard], end, last.id + 1 == manifest.blocks)
                    {
                        if let Ok(file) = readers[shard].lock() {
                            cache_calls += 1;
                            cache_requested += length;
                            if !release_consumed_spool(&file, offset, length) {
                                cache_failed += 1;
                            }
                        }
                        released[shard] = offset + length;
                    }
                }
            }
            for (id, acc, timing) in prepared {
                total.decode_match_seconds += timing.decode_match_seconds;
                total.assemble_seconds += timing.assemble_seconds;
                total.compression_seconds += timing.compression_seconds;
                total.blocks += 1;
                let before = Instant::now();
                tx.send((id, acc))?;
                total.queue_wait_seconds += before.elapsed().as_secs_f64();
            }
            if total.blocks % 100 == 0 {
                eprintln!(
                    "DROPLET_REPLAY_PROGRESS blocks={} total={} wall={:.3}",
                    total.blocks,
                    manifest.blocks,
                    start.elapsed().as_secs_f64()
                );
            }
        }
        Ok(())
    })();
    drop(tx);
    let _ = tx_completion.send(processing.is_ok());
    let written = writer
        .join()
        .map_err(|_| anyhow::anyhow!("replay writer panic"));
    let drained = drain.join().map_err(|_| anyhow::anyhow!("recycle panic"));
    written??;
    processing?;
    drained?;
    let report: serde_json::Value = serde_json::from_slice(&std::fs::read(&private_report)?)?;
    anyhow::ensure!(
        report["input_read_pairs"].as_u64() == Some(manifest.pairs),
        "final replay pair mismatch"
    );
    std::fs::rename(private_report, final_report)?;
    eprintln!(
        "DROPLET_REPLAY_PROFILE {}",
        serde_json::json!({
             "window":window, "wall_seconds":start.elapsed().as_secs_f64(),
            "spool_cache_advice_requested_bytes":cache_requested,
            "spool_cache_advice_calls":cache_calls, "spool_cache_advice_failed":cache_failed,
            "block_timings":total, "timing_note":"sums of block wall intervals overlap when window>1"
        })
    );
    Ok(())
}

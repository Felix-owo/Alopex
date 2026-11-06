use super::*;
use needletail::parse_fastx_file;
use std::time::Instant;

pub(super) fn write_fastq_record(
    output: &mut Vec<u8>,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
    trim_start: usize,
    umi: Option<&[u8]>,
    trim: bool,
    keep_comment: bool,
) -> Result<()> {
    output.push(b'@');
    if let Some(umi_seq) = umi {
        output.extend_from_slice(normalize_read_id(id)?);
        output.push(b':');
        output.extend_from_slice(umi_seq);
        if keep_comment {
            if let Some(comment_start) = id.iter().position(u8::is_ascii_whitespace) {
                output.extend_from_slice(&id[comment_start..]);
            }
        }
    } else {
        output.extend_from_slice(id);
    }
    output.push(b'\n');

    if trim {
        if trim_start < seq.len() {
            output.extend_from_slice(&seq[trim_start..]);
            output.extend_from_slice(b"\n+\n");
            output.extend_from_slice(&qual[trim_start..]);
        } else {
            output.extend_from_slice(b"N\n+\n#");
        }
    } else {
        output.extend_from_slice(seq);
        output.extend_from_slice(b"\n+\n");
        output.extend_from_slice(qual);
    }
    output.push(b'\n');
    Ok(())
}

pub(super) fn normalize_read_id(id: &[u8]) -> Result<&[u8]> {
    // needletail 的 id() 已移除 FASTQ marker；名称中的 @ 不能再删，
    // 否则不同分子会同名。首字符为空白时没有 ID，不能把备注提升成名称。
    let end = id
        .iter()
        .position(|base| base.is_ascii_whitespace())
        .unwrap_or(id.len());
    let mut normalized = &id[..end];
    if normalized.ends_with(b"/1") || normalized.ends_with(b"/2") {
        normalized = &normalized[..normalized.len() - 2];
    }
    if normalized.is_empty() {
        anyhow::bail!("FASTQ read identifier is empty after normalization");
    }
    Ok(normalized)
}

pub(super) fn validate_paired_read_ids<'a>(r1_id: &'a [u8], r2_id: &[u8]) -> Result<&'a [u8]> {
    let r1_normalized = normalize_read_id(r1_id)?;
    let r2_normalized = normalize_read_id(r2_id)?;
    if r1_normalized != r2_normalized {
        anyhow::bail!(
            "Paired-end read-ID mismatch: R1 '{}' does not match R2 '{}' after normalization",
            String::from_utf8_lossy(r1_id),
            String::from_utf8_lossy(r2_id)
        );
    }
    Ok(r1_normalized)
}

pub(super) fn validate_fastq_quality(seq: &[u8], qual: &[u8]) -> Result<()> {
    anyhow::ensure!(
        seq.len() == qual.len(),
        "FASTQ sequence/quality length mismatch"
    );
    anyhow::ensure!(
        qual.iter().all(|q| (33..=126).contains(q)),
        "Invalid FASTQ Phred+33 quality"
    );
    Ok(())
}

pub(super) fn read_fastq_pairs(
    r1_paths: &[PathBuf],
    r2_paths: &[PathBuf],
    chunk_size: usize,
    rx_input_recycle: &Receiver<InputChunk>,
    rx_acc_recycle: &Receiver<AccumulationBuffer>,
    tx_work: &Sender<WorkChunk>,
) -> Result<u64> {
    if r1_paths.is_empty() || r1_paths.len() != r2_paths.len() {
        anyhow::bail!(
            "FASTQ lane list mismatch: R1={} R2={}",
            r1_paths.len(),
            r2_paths.len()
        );
    }

    let mut current_chunk = rx_input_recycle
        .recv()
        .context("Input chunk pool closed before reading started")?;
    let mut current_acc = rx_acc_recycle
        .recv()
        .context("Accumulation buffer pool closed before reading started")?;
    let mut total = 0u64;
    let mut chunk_id = 0u64;
    let start_time = Instant::now();

    for (lane_index, (r1_path, r2_path)) in r1_paths.iter().zip(r2_paths.iter()).enumerate() {
        log::info!(
            "Reading FASTQ pair {}/{}: R1={} R2={}",
            lane_index + 1,
            r1_paths.len(),
            r1_path.display(),
            r2_path.display()
        );
        let mut reader1 = parse_fastx_file(r1_path)
            .with_context(|| format!("R1 open failed: {}", r1_path.display()))?;
        let mut reader2 = parse_fastx_file(r2_path)
            .with_context(|| format!("R2 open failed: {}", r2_path.display()))?;

        loop {
            match (reader1.next(), reader2.next()) {
                (Some(r1_result), Some(r2_result)) => {
                    let r1 = r1_result.with_context(|| {
                        format!("Error reading R1 record from {}", r1_path.display())
                    })?;
                    let r2 = r2_result.with_context(|| {
                        format!("Error reading R2 record from {}", r2_path.display())
                    })?;
                    validate_paired_read_ids(r1.id(), r2.id())?;
                    let r1_seq = r1.seq();
                    let is_usable = r1_seq.as_ref().len() >= MIN_R1_LEN;
                    let r2_seq = r2.seq();
                    let r1_qual = r1.qual().context("R1 record is missing quality scores")?;
                    let r2_qual = r2.qual().context("R2 record is missing quality scores")?;
                    validate_fastq_quality(r1_seq.as_ref(), r1_qual)?;
                    validate_fastq_quality(r2_seq.as_ref(), r2_qual)?;

                    current_chunk.push_pair(
                        r1.id(),
                        r1_seq.as_ref(),
                        r1_qual,
                        r2.id(),
                        r2_seq.as_ref(),
                        r2_qual,
                    );
                    if is_usable {
                        total += 1;
                    }

                    if current_chunk.count >= chunk_size {
                        tx_work
                            .send((chunk_id, current_chunk, current_acc))
                            .context("Worker channel closed while sending a full chunk")?;
                        chunk_id = chunk_id.checked_add(1).context("Input chunk id overflow")?;
                        current_chunk = rx_input_recycle
                            .recv()
                            .context("Input chunk recycle channel closed while reading")?;
                        current_acc = rx_acc_recycle
                            .recv()
                            .context("Accumulation buffer recycle channel closed while reading")?;
                    }

                    if is_usable && total % PROGRESS_EVERY_READS == 0 {
                        log::info!(
                            "Processed {:.1}M read pairs in {:.2}s",
                            total as f64 / 1e6,
                            start_time.elapsed().as_secs_f64()
                        );
                    }
                }
                (None, None) => break,
                (Some(extra), None) => {
                    let rec = extra.context("Error reading extra R1 record")?;
                    anyhow::bail!(
                        "Paired-end mismatch in FASTQ pair {}: R1 has extra record beyond R2; first extra id: {}",
                        lane_index + 1,
                        String::from_utf8_lossy(rec.id())
                    );
                }
                (None, Some(extra)) => {
                    let rec = extra.context("Error reading extra R2 record")?;
                    anyhow::bail!(
                        "Paired-end mismatch in FASTQ pair {}: R2 has extra record beyond R1; first extra id: {}",
                        lane_index + 1,
                        String::from_utf8_lossy(rec.id())
                    );
                }
            }
        }
    }

    if current_chunk.count > 0 {
        tx_work
            .send((chunk_id, current_chunk, current_acc))
            .context("Worker channel closed while sending the final chunk")?;
    }

    log::info!(
        "Finished reading {} usable read pairs across {} FASTQ pair(s)",
        total,
        r1_paths.len()
    );
    Ok(total)
}

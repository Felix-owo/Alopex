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
) -> Result<()> {
    output.push(b'@');
    if let Some(umi_seq) = umi {
        let id = parse_read_id(id)?;
        output.extend_from_slice(id.name);
        if let Some(barcode) = id.barcode {
            output.push(b':');
            output.extend_from_slice(barcode);
        }
        output.push(b':');
        output.extend_from_slice(umi_seq);
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

struct ReadId<'a> {
    name: &'a [u8],
    mate: Option<u8>,
    barcode: Option<&'a [u8]>,
}

pub(super) fn is_read_id_whitespace(&byte: &u8) -> bool {
    matches!(byte, b' ' | b'\t'..=b'\r')
}

fn illumina_mate(field: &[u8]) -> Option<u8> {
    let [side @ (b'1' | b'2'), b':', b'Y' | b'N', b':', rest @ ..] = field else {
        return None;
    };
    let colon = rest.iter().position(|&b| b == b':')?;
    let control = &rest[..colon];
    let index = &rest[colon + 1..];
    if control.is_empty() || !control.iter().all(u8::is_ascii_digit) {
        return None;
    }
    // 末字段可以是样本号、单/双 index，或留空；只取读段端，不将 index 当成配对键。
    let sample_number = index.iter().all(u8::is_ascii_digit);
    let is_index = |s: &[u8]| !s.is_empty() && s.iter().all(|b| b"ACGTN".contains(b));
    (sample_number || index.splitn(2, |&b| b == b'+').all(is_index)).then_some(*side)
}

fn parse_read_id(id: &[u8]) -> Result<ReadId<'_>> {
    // needletail 的 id() 已移除 FASTQ marker；名称中的 @ 不能再删，
    // 否则不同分子会同名。首字符为空白时没有 ID，不能把备注提升成名称。
    let mut fields = id.split(is_read_id_whitespace);
    let mut name = fields.next().unwrap_or_default();
    let mut mate = None;
    let mut barcode = None;
    if name.ends_with(b"/1") || name.ends_with(b"/2") {
        mate = name.last().copied();
        name = &name[..name.len() - 2];
    } else if let Some(offset) = name.len().checked_sub(3 + 17) {
        // SeekGene 的 mate 标记在 CB17 注释之前；注释仍是配对身份的一部分。
        let (stem, suffix) = name.split_at(offset);
        if let [b'/', side @ (b'1' | b'2'), b':', cb @ ..] = suffix {
            if cb.iter().all(|base| b"ACGTN".contains(base)) {
                name = stem;
                mate = Some(*side);
                barcode = Some(cb);
            }
        }
    }
    if name.is_empty() {
        anyhow::bail!("FASTQ read identifier is empty after normalization");
    }
    if let Some(side) = fields.find(|field| !field.is_empty()).and_then(illumina_mate) {
        anyhow::ensure!(
            mate.is_none_or(|suffix| suffix == side),
            "Conflicting FASTQ mate markers in '{}'", String::from_utf8_lossy(id)
        );
        mate = Some(side);
    }
    Ok(ReadId { name, mate, barcode })
}

pub(super) fn validate_paired_read_ids(r1_id: &[u8], r2_id: &[u8]) -> Result<()> {
    let r1 = parse_read_id(r1_id)?;
    let r2 = parse_read_id(r2_id)?;
    if r1.name != r2.name || r1.barcode != r2.barcode
        || r1.mate == Some(b'2') || r2.mate == Some(b'1') {
        anyhow::bail!(
            "Paired-end read-ID mismatch: R1 '{}' does not match R2 '{}' after normalization",
            String::from_utf8_lossy(r1_id),
            String::from_utf8_lossy(r2_id)
        );
    }
    Ok(())
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

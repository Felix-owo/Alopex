mod barcode_map {

    use super::*;

    pub(super) fn load_barcode_map(path: &Path, mode: DemuxMode) -> Result<Vec<CellRow>> {
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .delimiter(b',')
            .from_path(path)
            .with_context(|| format!("Cannot read Barcode_Map.csv: {}", path.display()))?;

        let headers = reader
            .headers()
            .context("Cannot read Barcode_Map.csv header")?
            .clone();
        let mut seen_headers = AHashSet::new();
        for header in &headers {
            let key = header
                .trim_start_matches('\u{feff}')
                .trim()
                .to_ascii_lowercase();
            if key.is_empty() || !seen_headers.insert(key) {
                anyhow::bail!(
                    "Barcode_Map.csv has a blank or duplicate column: {:?}",
                    header
                );
            }
        }
        let idx_dna = header_index(&headers, "DNA_Barcode")?;
        let idx_rna = header_index_optional(&headers, "RNA_Barcode");
        if mode == DemuxMode::DnaRna && idx_rna.is_none() {
            anyhow::bail!("Barcode_Map.csv missing required column RNA_Barcode for dna-rna mode");
        }
        let idx_plate = header_index_optional(&headers, "PlateID");
        if mode != DemuxMode::DnaOnlyDroplet && idx_plate.is_none() {
            anyhow::bail!("Barcode_Map.csv missing required column PlateID");
        }
        let idx_order = header_index(&headers, "Cell_Order")?;

        let mut rows = Vec::new();
        for (line_idx, record) in reader.records().enumerate() {
            let record = record.with_context(|| {
                format!(
                    "Cannot parse Barcode_Map.csv record at data line {}",
                    line_idx + 1
                )
            })?;
            let dna = record.get(idx_dna).unwrap_or("").trim().to_string();
            let rna = idx_rna
                .and_then(|idx| record.get(idx))
                .unwrap_or("")
                .trim()
                .to_string();
            let plate = idx_plate
                .and_then(|i| record.get(i))
                .unwrap_or("")
                .trim()
                .to_string();
            let order = record.get(idx_order).unwrap_or("").trim().to_string();

            if dna.is_empty() && rna.is_empty() && plate.is_empty() && order.is_empty() {
                continue;
            }
            if dna.is_empty()
                || (mode != DemuxMode::DnaOnlyDroplet && plate.is_empty())
                || order.is_empty()
                || (mode == DemuxMode::DnaRna && rna.is_empty())
            {
                anyhow::bail!(
                "Barcode_Map.csv row {} has empty required field(s): DNA_Barcode='{}', RNA_Barcode='{}', PlateID='{}', Cell_Order='{}'",
                line_idx + 2,
                dna,
                rna,
                plate,
                order
            );
            }
            if mode == DemuxMode::DnaOnlyDroplet {
                if dna.len() != 17
                    || !dna.bytes().all(|b| matches!(b, b'A' | b'G' | b'T'))
                    || !plate.is_empty()
                    || !rna.is_empty()
                    || !droplet::design_barcodes()?.contains(dna.as_bytes())
                {
                    anyhow::bail!("Droplet called cells require fixed-design 17bp DNA_Barcode and empty plate/RNA identities");
                }
            } else {
                validate_barcode_alphabet(&dna, "DNA_Barcode", line_idx + 2)?;
            }
            if matches!(plate.as_str(), "." | "..")
                || !plate
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-'))
            {
                anyhow::bail!(
                "Barcode_Map.csv row {} PlateID '{}' contains unsupported path characters; use only letters, digits, '.', '_' or '-'",
                line_idx + 2,
                plate
            );
            }
            if !rna.is_empty() {
                validate_barcode_alphabet(&rna, "RNA_Barcode", line_idx + 2)?;
            }
            rows.push(CellRow {
                dna_barcode: dna,
                rna_barcode: rna,
                plate_id: plate,
                cell_order: order,
            });
        }

        if rows.is_empty() {
            anyhow::bail!(
                "Barcode_Map.csv contains no usable rows: {}",
                path.display()
            );
        }
        validate_barcode_rows(&rows, mode)?;
        Ok(rows)
    }

    pub(super) fn header_index(headers: &csv::StringRecord, name: &str) -> Result<usize> {
        header_index_optional(headers, name)
            .with_context(|| format!("Barcode_Map.csv missing required column {}", name))
    }

    pub(super) fn header_index_optional(headers: &csv::StringRecord, name: &str) -> Option<usize> {
        headers.iter().position(|header| {
            header
                .trim_start_matches('\u{feff}')
                .trim()
                .eq_ignore_ascii_case(name)
        })
    }

    pub(super) fn validate_barcode_alphabet(
        barcode: &str,
        column: &str,
        row_number: usize,
    ) -> Result<()> {
        if let Some(invalid) = barcode
            .bytes()
            .find(|base| !matches!(*base, b'A' | b'C' | b'G' | b'T'))
        {
            anyhow::bail!(
            "Barcode_Map.csv row {} column {} contains invalid base '{}'. Protected barcodes must use raw uppercase A/C/G/T",
            row_number,
            column,
            invalid as char
        );
        }
        if !matches!(barcode.len(), 8 | 10) {
            anyhow::bail!(
                "Barcode_Map.csv row {} column {} must be 8 or 10 bases, found {} for '{}'",
                row_number,
                column,
                barcode.len(),
                barcode
            );
        }
        Ok(())
    }

    pub(super) fn validate_barcode_rows(rows: &[CellRow], mode: DemuxMode) -> Result<()> {
        let mut dna_map: AHashMap<String, (String, String)> = AHashMap::new();
        let mut rna_map: AHashMap<String, (String, String)> = AHashMap::new();
        let mut plate_map: AHashMap<String, String> = AHashMap::new();
        let mut order_map: AHashMap<String, String> = AHashMap::new();

        for row in rows {
            validate_unique_barcode(
                &mut dna_map,
                "DNA_Barcode",
                &row.dna_barcode,
                &row.plate_id,
                &row.cell_order,
            )?;
            if mode == DemuxMode::DnaRna {
                validate_unique_barcode(
                    &mut rna_map,
                    "RNA_Barcode",
                    &row.rna_barcode,
                    &row.plate_id,
                    &row.cell_order,
                )?;
            }

            if let Some(existing) = (!row.plate_id.is_empty())
                .then(|| plate_map.insert(row.plate_id.clone(), row.dna_barcode.clone()))
                .flatten()
            {
                anyhow::bail!(
                    "Barcode_Map.csv PlateID '{}' is duplicated for DNA barcodes '{}' and '{}'",
                    row.plate_id,
                    existing,
                    row.dna_barcode
                );
            }
            if let Some(existing) =
                order_map.insert(row.cell_order.clone(), row.dna_barcode.clone())
            {
                anyhow::bail!(
                    "Barcode_Map.csv Cell_Order '{}' is duplicated for DNA barcodes '{}' and '{}'",
                    row.cell_order,
                    existing,
                    row.dna_barcode
                );
            }
        }
        if mode != DemuxMode::DnaOnlyDroplet {
            audit_correction_neighborhoods(dna_map.keys(), "DNA_Barcode")?;
        }
        if mode == DemuxMode::DnaRna {
            audit_correction_neighborhoods(rna_map.keys(), "RNA_Barcode")?;
        }
        Ok(())
    }

    fn audit_correction_neighborhoods<'a>(
        barcodes: impl Iterator<Item = &'a String>,
        label: &str,
    ) -> Result<()> {
        let values = barcodes.collect::<Vec<_>>();
        for (i, left) in values.iter().enumerate() {
            for right in values.iter().skip(i + 1) {
                if left.len() != right.len() {
                    continue;
                }
                let distance = left
                    .bytes()
                    .zip(right.bytes())
                    .filter(|(a, b)| a != b)
                    .count();
                if distance <= 2 {
                    anyhow::bail!(
                    "{} one-substitution correction neighborhoods overlap: '{}' and '{}' have Hamming distance {}",
                    label,
                    left,
                    right,
                    distance
                );
                }
            }
        }
        Ok(())
    }

    pub(super) fn validate_unique_barcode(
        map: &mut AHashMap<String, (String, String)>,
        column: &str,
        barcode: &str,
        plate_id: &str,
        cell_order: &str,
    ) -> Result<()> {
        if let Some(entry) = map.insert(
            barcode.to_string(),
            (plate_id.to_string(), cell_order.to_string()),
        ) {
            anyhow::bail!(
                "Barcode_Map.csv duplicate {} '{}': rows map to ({}, {}) and ({}, {})",
                column,
                barcode,
                entry.0,
                entry.1,
                plate_id,
                cell_order
            );
        }
        Ok(())
    }

    pub(super) fn resolve_barcode_lengths(
        rows: &[CellRow],
        mode: DemuxMode,
    ) -> Result<(Vec<usize>, usize)> {
        let dna_lengths = infer_barcode_lengths(rows, "DNA", |row| row.dna_barcode.as_str())?;
        let rna_len = if mode == DemuxMode::DnaRna {
            infer_uniform_barcode_len(rows, "RNA", |row| row.rna_barcode.as_str())?
        } else {
            0
        };
        Ok((dna_lengths, rna_len))
    }

    pub(super) fn infer_barcode_lengths<F>(
        rows: &[CellRow],
        label: &str,
        get_barcode: F,
    ) -> Result<Vec<usize>>
    where
        F: Fn(&CellRow) -> &str,
    {
        let mut lengths = rows
            .iter()
            .map(|row| get_barcode(row).len())
            .collect::<Vec<_>>();
        lengths.sort_unstable_by(|a, b| b.cmp(a));
        lengths.dedup();
        match lengths.as_slice() {
            [] => anyhow::bail!("Barcode_Map.csv has no {} barcode rows", label),
            [0] => anyhow::bail!("{} barcode length cannot be 0", label),
            _ if lengths.contains(&0) => anyhow::bail!("{} barcode length cannot be 0", label),
            _ => Ok(lengths),
        }
    }

    pub(super) fn infer_uniform_barcode_len<F>(
        rows: &[CellRow],
        label: &str,
        get_barcode: F,
    ) -> Result<usize>
    where
        F: Fn(&CellRow) -> &str,
    {
        let lengths = infer_barcode_lengths(rows, label, get_barcode)?;
        match lengths.as_slice() {
        [len] => Ok(*len),
        _ => anyhow::bail!(
            "{} barcode lengths are not uniform: {:?}. Mixed RNA lengths are not supported because TSO, UMI, and trim positions become ambiguous.",
            label,
            lengths
        ),
    }
    }

    pub(super) fn validate_effective_barcode_prefixes(
        rows: &[CellRow],
        dna_lengths: &[usize],
        rna_len: usize,
        mode: DemuxMode,
    ) -> Result<()> {
        let mixed_dna_lengths = dna_lengths.len() > 1;
        for &dna_len in dna_lengths {
            validate_effective_prefixes_for_modality(
                rows,
                Modality::Dna,
                dna_len,
                mixed_dna_lengths.then_some(dna_len),
                |row| row.dna_barcode.as_str(),
            )?;
        }
        if mode == DemuxMode::DnaRna {
            validate_effective_prefixes_for_modality(rows, Modality::Rna, rna_len, None, |row| {
                row.rna_barcode.as_str()
            })?;
        }
        Ok(())
    }

    pub(super) fn validate_effective_prefixes_for_modality<F>(
        rows: &[CellRow],
        modality: Modality,
        barcode_len: usize,
        actual_len_filter: Option<usize>,
        get_barcode: F,
    ) -> Result<()>
    where
        F: Fn(&CellRow) -> &str,
    {
        let mut seen: AHashMap<String, String> = AHashMap::new();
        for row in rows {
            let barcode = get_barcode(row);
            if actual_len_filter.is_some_and(|actual_len| barcode.len() != actual_len) {
                continue;
            }
            let raw_prefix = barcode[..barcode_len].to_string();
            if let Some(existing) = seen.insert(raw_prefix.clone(), barcode.to_string()) {
                if existing != barcode {
                    anyhow::bail!(
                    "{} raw barcode prefix collision at length {}: '{}' and '{}' both use effective prefix '{}'. Increase barcode length or fix Barcode_Map.csv.",
                    modality.as_str(),
                    barcode_len,
                    existing,
                    barcode,
                    raw_prefix
                );
                }
            }
        }
        Ok(())
    }
}
mod cli {

    use clap::{Parser, ValueEnum};
    use std::path::PathBuf;

    #[derive(Parser, Debug, Clone)]
    #[command(author, version, about = "DNA/RNA FASTQ 高通量解复用器")]
    pub(crate) struct Args {
        #[arg(
            long,
            required = true,
            help = "R1 FASTQ；可按 lane/chunk 顺序重复传入。"
        )]
        pub(crate) r1: Vec<PathBuf>,

        #[arg(
            long,
            required_unless_present_any = ["count_only", "evidence_pairs"],
            help = "R2 FASTQ；数量和顺序必须与 --r1 一一对应。"
        )]
        pub(crate) r2: Vec<PathBuf>,

        #[arg(long, requires = "counts_output")]
        pub(crate) count_only: bool,
        #[arg(long, requires = "count_only")]
        pub(crate) counts_output: Option<PathBuf>,
        #[arg(long, requires = "count_only")]
        pub(crate) count_metrics: Option<PathBuf>,
        #[arg(long, requires = "evidence_output", conflicts_with = "count_only")]
        pub(crate) evidence_pairs: Option<PathBuf>,
        #[arg(long, requires = "evidence_pairs")]
        pub(crate) evidence_output: Option<PathBuf>,

        #[arg(
            long = "barcode-map",
            help = "Barcode_Map.csv；必需 DNA_Barcode、PlateID、Cell_Order，dna-rna 模式还需 RNA_Barcode。"
        )]
        pub(crate) barcode_map: Option<PathBuf>,

        #[arg(long = "dna-out-dir", help = "DNA demux FASTQ 输出根目录。")]
        pub(crate) dna_out_dir: Option<PathBuf>,

        #[arg(long = "rna-out-dir", help = "RNA demux FASTQ 输出根目录。")]
        pub(crate) rna_out_dir: Option<PathBuf>,

        #[arg(long = "json-report", help = "JSON report 输出路径。")]
        pub(crate) json_report: Option<PathBuf>,

        #[arg(long = "sample-name", allow_hyphen_values = true, default_value = "")]
        pub(crate) sample_name: String,

        #[arg(
            long = "min-matched-read-pairs",
            default_value_t = 1,
            help = "保留 barcode bucket 所需的最少 matched read pairs。"
        )]
        pub(crate) min_matched_read_pairs: u64,

        #[arg(long, default_value_t = 8)]
        pub(crate) threads: usize,

        #[arg(
            long = "dna-w-spacer-len",
            default_value_t = 0,
            help = "DNA barcode 与 ME motif 之间固定 spacer 的碱基数。"
        )]
        pub(crate) dna_w_spacer_len: usize,

        #[arg(
        long = "mode",
        value_enum,
        default_value_t = DemuxMode::DnaOnly,
        help = "dna-only 用于 Cabernet；dna-only-taps 使用普通四碱基 ME；dna-rna 用于 SRD DNA-ME/TSO-RNA 混合布局。"
    )]
        pub(crate) mode: DemuxMode,
    }

    #[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
    pub(crate) enum DemuxMode {
        #[value(name = "dna-only")]
        DnaOnly,
        #[value(name = "dna-only-taps")]
        DnaOnlyTaps,
        #[value(name = "dna-only-droplet")]
        DnaOnlyDroplet,
        #[value(name = "dna-rna")]
        DnaRna,
    }
}
mod fastq_io {

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
            output.extend_from_slice(normalize_read_id(id)?);
            output.push(b':');
            output.extend_from_slice(umi_seq);
            let header = id.trim_ascii_start();
            if let Some(comment_start) = header.iter().position(u8::is_ascii_whitespace) {
                output.extend_from_slice(&header[comment_start..]);
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
        let mut normalized = id;
        while normalized.first().is_some_and(u8::is_ascii_whitespace) {
            normalized = &normalized[1..];
        }
        if normalized.first() == Some(&b'@') {
            normalized = &normalized[1..];
        }
        let end = normalized
            .iter()
            .position(|base| base.is_ascii_whitespace())
            .unwrap_or(normalized.len());
        normalized = &normalized[..end];
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
                            chunk_id =
                                chunk_id.checked_add(1).context("Input chunk id overflow")?;
                            current_chunk = rx_input_recycle
                                .recv()
                                .context("Input chunk recycle channel closed while reading")?;
                            current_acc = rx_acc_recycle.recv().context(
                                "Accumulation buffer recycle channel closed while reading",
                            )?;
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
}
mod matcher {

    use super::*;

    pub(super) fn process_read<'a>(
        r1: &'a [u8],
        dna_indices: &'a [BarcodeIndex],
        rna_idx: &'a BarcodeIndex,
        dna_lut: &'static [u8; 256],
        rna_lut: &'static [u8; 256],
        dna_w_spacer_len: usize,
        demux_mode: DemuxMode,
    ) -> ProcResult<'a> {
        if demux_mode == DemuxMode::DnaOnlyDroplet {
            return droplet::match_read(r1, &dna_indices[0]);
        }
        let mut dna_match: Option<ProcResult<'a>> = None;
        let mut dna_rejection: Option<ProcResult<'a>> = None;
        for dna_idx in dna_indices {
            if let Some(result) = match_dna_read(r1, dna_idx, dna_lut, dna_w_spacer_len) {
                if result.is_valid {
                    if dna_match.is_some() {
                        return ProcResult::cross_layout_ambiguous();
                    }
                    dna_match = Some(result);
                } else {
                    retain_best_rejection(&mut dna_rejection, result);
                }
            }
        }
        if let Some(result) = dna_match {
            return result;
        }

        if demux_mode != DemuxMode::DnaRna {
            return dna_rejection.unwrap_or_else(ProcResult::unmatched);
        }

        let rna_rejection = match_rna_read(r1, rna_idx, rna_lut);
        if let Some(result) = rna_rejection {
            if result.is_valid {
                if dna_rejection.is_some() {
                    return ProcResult::cross_layout_ambiguous();
                }
                return result;
            }
            if result.fate == TerminalFate::Ambiguous {
                return result;
            }
            if dna_rejection
                .as_ref()
                .is_some_and(|candidate| candidate.fate == TerminalFate::Ambiguous)
            {
                return dna_rejection.unwrap();
            }
            return result;
        }
        dna_rejection.unwrap_or_else(ProcResult::unmatched)
    }

    pub(super) fn retain_best_rejection<'a>(
        slot: &mut Option<ProcResult<'a>>,
        candidate: ProcResult<'a>,
    ) {
        let should_replace = slot.is_none()
            || (candidate.fate == TerminalFate::Ambiguous
                && slot
                    .as_ref()
                    .is_some_and(|current| current.fate != TerminalFate::Ambiguous));
        if should_replace {
            *slot = Some(candidate);
        }
    }

    pub(super) fn match_rna_read<'a>(
        r1: &'a [u8],
        rna_idx: &'a BarcodeIndex,
        rna_lut: &'static [u8; 256],
    ) -> Option<ProcResult<'a>> {
        let rna_len = rna_idx.len();
        let mut rejection = None;
        let search_window = if r1.len() > 60 { &r1[..60] } else { r1 };

        for offset in 0..=1 {
            if offset + TSO_SEQ.len() <= r1.len() {
                let anchor_distance =
                    lut_distance(TSO_SEQ, &r1[offset..offset + TSO_SEQ.len()], rna_lut);
                if anchor_distance > 2 {
                    continue;
                }
                let barcode_start = TSO_SEQ.len() + offset;
                let umi_start = barcode_start + rna_len;
                if umi_start <= r1.len() {
                    let observed = &r1[barcode_start..umi_start];
                    let barcode_match = rna_idx.match_simple(observed);
                    let umi_end = umi_start + LEN_UMI;
                    let payload_start = umi_end + LEN_SPACER;
                    let raw_umi = if umi_end <= r1.len() {
                        Some(&r1[umi_start..umi_end])
                    } else {
                        None
                    };
                    if payload_start < r1.len() {
                        if let (BarcodeMatch::Unique { id }, Some(umi)) = (barcode_match, raw_umi) {
                            return Some(ProcResult::assigned(
                                Modality::Rna,
                                id,
                                payload_start,
                                Some(umi),
                            ));
                        }
                    }
                    retain_best_rejection(
                        &mut rejection,
                        rejected_barcode_attempt(barcode_match, raw_umi),
                    );
                }
            }
        }

        if let Some(alignment) = align_fallback_detail(TSO_SEQ, search_window, rna_lut) {
            let umi_start = alignment.end + rna_len;
            if umi_start <= r1.len() {
                let observed = &r1[alignment.end..umi_start];
                let barcode_match = rna_idx.match_simple(observed);
                let umi_end = umi_start + LEN_UMI;
                let payload_start = umi_end + LEN_SPACER;
                let raw_umi = if umi_end <= r1.len() {
                    Some(&r1[umi_start..umi_end])
                } else {
                    None
                };
                if payload_start < r1.len() {
                    if let (BarcodeMatch::Unique { id }, Some(umi)) = (barcode_match, raw_umi) {
                        return Some(ProcResult::assigned(
                            Modality::Rna,
                            id,
                            payload_start,
                            Some(umi),
                        ));
                    }
                }
                retain_best_rejection(
                    &mut rejection,
                    rejected_barcode_attempt(barcode_match, raw_umi),
                );
            }
        }
        rejection
    }

    pub(super) fn match_dna_read<'a>(
        r1: &'a [u8],
        dna_idx: &'a BarcodeIndex,
        dna_lut: &'static [u8; 256],
        dna_w_spacer_len: usize,
    ) -> Option<ProcResult<'a>> {
        let dna_len = dna_idx.len();
        let mut rejection = None;

        let dna_me_pos = dna_len + dna_w_spacer_len;
        let direct_attempts = [
            (dna_me_pos, 0usize, dna_len),
            (
                dna_len.saturating_sub(1) + dna_w_spacer_len,
                0usize,
                dna_len.saturating_sub(1),
            ),
            (dna_len + 1 + dna_w_spacer_len, 1usize, dna_len + 1),
        ];
        for (anchor_offset, barcode_start, barcode_end) in direct_attempts {
            if barcode_end <= barcode_start
                || barcode_end > r1.len()
                || anchor_offset + ME_SEQ.len() > r1.len()
            {
                continue;
            }
            let anchor_distance = lut_distance(
                ME_SEQ,
                &r1[anchor_offset..anchor_offset + ME_SEQ.len()],
                dna_lut,
            );
            if anchor_distance > 2 {
                continue;
            }
            let observed = &r1[barcode_start..barcode_end];
            let barcode_match = dna_idx.match_cascade(observed);
            let payload_start = anchor_offset + ME_SEQ.len() + LEN_GAP_AFTER_ME;
            if payload_start < r1.len() {
                if let BarcodeMatch::Unique { id } = barcode_match {
                    return Some(ProcResult::assigned(Modality::Dna, id, payload_start, None));
                }
            }
            retain_best_rejection(
                &mut rejection,
                rejected_barcode_attempt(barcode_match, None),
            );
        }

        let search_window = if r1.len() > 60 { &r1[..60] } else { r1 };
        if let Some(alignment) = align_fallback_detail(ME_SEQ, search_window, dna_lut) {
            let pos = alignment.start;
            let expected_me_pos = dna_len + dna_w_spacer_len;
            let lower = expected_me_pos.saturating_sub(1);
            let upper = expected_me_pos + 1;
            if pos >= lower && pos <= upper && pos >= dna_w_spacer_len {
                let barcode_end = pos - dna_w_spacer_len;
                let barcode_start = barcode_end.saturating_sub(dna_len);
                if barcode_start <= barcode_end && barcode_end <= r1.len() {
                    let observed = &r1[barcode_start..barcode_end];
                    let barcode_match = dna_idx.match_cascade(observed);
                    let payload_start = pos + ME_SEQ.len() + LEN_GAP_AFTER_ME;
                    if payload_start < r1.len() {
                        if let BarcodeMatch::Unique { id } = barcode_match {
                            return Some(ProcResult::assigned(
                                Modality::Dna,
                                id,
                                payload_start,
                                None,
                            ));
                        }
                    }
                    retain_best_rejection(
                        &mut rejection,
                        rejected_barcode_attempt(barcode_match, None),
                    );
                }
            }
        }

        rejection
    }

    pub(super) fn rejected_barcode_attempt<'a>(
        barcode_match: BarcodeMatch<'_>,
        umi: Option<&'a [u8]>,
    ) -> ProcResult<'a> {
        let fate = match barcode_match {
            BarcodeMatch::Ambiguous => TerminalFate::Ambiguous,
            BarcodeMatch::Unique { .. } | BarcodeMatch::None => TerminalFate::Unmatched,
        };
        ProcResult::rejected(fate, umi)
    }

    pub(super) fn lut_distance(query: &[u8], target: &[u8], lut: &[u8; 256]) -> usize {
        query
            .iter()
            .zip(target.iter())
            .filter(|(query_base, target_base)| {
                lut[**query_base as usize] != lut[**target_base as usize]
            })
            .count()
            + query.len().abs_diff(target.len())
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub(super) struct AnchorAlignment {
        pub(super) start: usize,
        pub(super) end: usize,
    }

    pub(super) fn align_fallback_detail(
        query: &[u8],
        target: &[u8],
        lut: &[u8; 256],
    ) -> Option<AnchorAlignment> {
        use triple_accel::levenshtein::levenshtein_search;
        let query_converted: Vec<u8> = query.iter().map(|&base| lut[base as usize]).collect();
        let target_converted: Vec<u8> = target.iter().map(|&base| lut[base as usize]).collect();
        levenshtein_search(&query_converted, &target_converted)
            .filter(|match_result| match_result.k <= 3)
            .min_by_key(|match_result| match_result.k)
            .map(|match_result| AnchorAlignment {
                start: match_result.start,
                end: match_result.end,
            })
    }
}
mod output {

    use super::*;
    use std::collections::VecDeque;
    use std::fs::{File, OpenOptions};
    use std::io::{BufWriter, Write};
    use std::process::Command;

    pub(super) fn compress_accumulation(acc: &mut AccumulationBuffer) -> Result<()> {
        let mut compressor = Compressor::new(CompressionLvl::fastest());
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

    pub(super) fn compress_plain_buffer(
        compressor: &mut Compressor,
        input: &[u8],
        output: &mut Vec<u8>,
        label: &str,
    ) -> Result<()> {
        let bound = compressor.gzip_compress_bound(input.len());
        if output.len() < bound {
            output.resize(bound, 0);
        }
        let compressed_len = compressor
            .gzip_compress(input, output)
            .with_context(|| format!("gzip compression failed for {label}"))?;
        output.truncate(compressed_len);
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
    ) -> JoinHandle<Result<u64>> {
        std::thread::spawn(move || writer_loop(rx_write, tx_acc_recycle, paths, cell_rows))
    }

    struct WriterCache {
        writers: AHashMap<BucketKey, (BufWriter<File>, BufWriter<File>)>,
        order: VecDeque<BucketKey>,
        max_pairs: usize,
        all_buckets_fit: bool,
    }

    const MAX_OPEN_BUCKETS: usize = 128;
    const MAX_RETAINED_DROPLET_BUCKETS: usize = 8192;
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
        let requested = if mode == DemuxMode::DnaOnlyDroplet
            && buckets <= MAX_RETAINED_DROPLET_BUCKETS
            && soft_limit.is_some()
        {
            buckets.max(MAX_OPEN_BUCKETS)
        } else {
            MAX_OPEN_BUCKETS
        };
        let available = safe_writer_pair_limit(requested, soft_limit);
        if available < buckets {
            available.min(MAX_OPEN_BUCKETS)
        } else {
            available
        }
    }

    #[test]
    fn droplet_writer_budget_retains_all_or_uses_bounded_fallback() {
        let mode = DemuxMode::DnaOnlyDroplet;
        assert_eq!(writer_pair_budget(mode, 6000, Some(12064)), 6000);
        assert_eq!(writer_pair_budget(mode, 6000, Some(12063)), 128);
        assert_eq!(writer_pair_budget(mode, 6000, Some(256)), 96);
        assert_eq!(writer_pair_budget(mode, 6000, None), 128);
        assert_eq!(writer_pair_budget(mode, 9000, Some(100000)), 128);
        assert_eq!(
            writer_pair_budget(DemuxMode::DnaRna, 6000, Some(100000)),
            128
        );
    }

    impl WriterCache {
        fn new(mode: DemuxMode, buckets: usize) -> Self {
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
                order: VecDeque::new(),
                max_pairs: effective_max_pairs,
                all_buckets_fit: mode == DemuxMode::DnaOnlyDroplet
                    && buckets <= effective_max_pairs,
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
        ) -> Result<&mut (BufWriter<File>, BufWriter<File>)> {
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
                self.writers
                    .insert(key.clone(), create_fastq_writers(key, paths)?);
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
    ) -> Result<u64> {
        let mut writers = WriterCache::new(
            paths.mode,
            paths.dna_cell_ids.len() + paths.rna_cell_ids.len(),
        );
        let mut stats: AHashMap<BucketKey, u64> = AHashMap::new();
        let mut matched_reads = 0u64;
        let mut fate_counts = FateCounts::default();

        let mut pending = BTreeMap::new();
        let mut next_chunk_id = 0u64;
        for (chunk_id, buffer) in rx_write {
            for mut buffer in
                insert_completed_chunk(&mut pending, &mut next_chunk_id, chunk_id, buffer)?
            {
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

        writers.flush_all()?;

        let fate_matched = fate_counts.matched_read_pairs();
        if fate_matched != matched_reads {
            anyhow::bail!(
            "Internal demux accounting mismatch: bucket matched_reads={} but fate matched_reads={}",
            matched_reads,
            fate_matched
        );
        }
        write_report(&paths, &cell_rows, &stats, fate_counts)?;
        Ok(matched_reads)
    }

    fn create_fastq_writers(
        key: &BucketKey,
        paths: &OutputPaths,
    ) -> Result<(BufWriter<File>, BufWriter<File>)> {
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
        Ok((BufWriter::new(r1), BufWriter::new(r2)))
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
                .with_context(|| {
                    format!("No canonical Sample ID for DNA barcode {}", row.dna_barcode)
                })?
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
        let file = File::create(&paths.json_report).with_context(|| {
            format!("Cannot create JSON report: {}", paths.json_report.display())
        })?;
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

    pub(super) fn prepare_output_roots(config: &Config) -> Result<()> {
        let dna_root = materialize_output_root(&config.dna_out_dir, "--dna-out-dir")?;
        let mut roots = vec![("--dna-out-dir", dna_root.clone())];
        if config.mode == DemuxMode::DnaRna {
            let rna_root = materialize_output_root(&config.rna_out_dir, "--rna-out-dir")?;
            if dna_root.starts_with(&rna_root) || rna_root.starts_with(&dna_root) {
                anyhow::bail!(
                    "--dna-out-dir and --rna-out-dir must be distinct, non-nested run roots"
                );
            }
            roots.push(("--rna-out-dir", rna_root));
        }
        for (label, root) in &roots {
            for (input_label, input_path) in config
                .r1
                .iter()
                .map(|path| ("--r1", path))
                .chain(config.r2.iter().map(|path| ("--r2", path)))
                .chain(std::iter::once(("--barcode-map", &config.barcode_map)))
            {
                let input = input_path.canonicalize().with_context(|| {
                    format!("Cannot resolve {input_label}: {}", input_path.display())
                })?;
                if input.starts_with(root) {
                    anyhow::bail!(
                        "Refusing to reset {label} because it contains {input_label}: {}",
                        input.display()
                    );
                }
            }
            if root.exists() {
                let metadata = std::fs::symlink_metadata(root)
                    .with_context(|| format!("Cannot inspect {label}: {}", root.display()))?;
                if metadata.file_type().is_symlink() || !metadata.is_dir() {
                    anyhow::bail!("{label} must be a real directory: {}", root.display());
                }
                std::fs::remove_dir_all(root).with_context(|| {
                    format!(
                        "Cannot reset exclusive run output root {label}: {}",
                        root.display()
                    )
                })?;
            }
            std::fs::create_dir_all(root)
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
}

mod droplet {
    use super::*;
    use std::io::{BufRead, BufReader, Write};

    static DESIGN_BARCODES: OnceLock<AHashSet<Vec<u8>>> = OnceLock::new();
    const DESIGN_BYTES: &[u8] =
        include_bytes!("../../../resources/droplet_ME5_U3CB_methylation.txt.gz");

    pub(super) fn design_barcodes() -> Result<&'static AHashSet<Vec<u8>>> {
        if let Some(design) = DESIGN_BARCODES.get() {
            return Ok(design);
        }
        let mut design = AHashSet::with_capacity(829_440);
        for line in BufReader::new(flate2::read::GzDecoder::new(DESIGN_BYTES)).lines() {
            let barcode = line?.into_bytes();
            anyhow::ensure!(
                barcode.len() == 17 && barcode.iter().all(|b| b"AGT".contains(b)),
                "Invalid DD-MET5 design barcode"
            );
            anyhow::ensure!(design.insert(barcode), "Duplicate DD-MET5 design barcode");
        }
        anyhow::ensure!(
            design.len() == 829_440,
            "Incomplete DD-MET5 design barcode resource"
        );
        let _ = DESIGN_BARCODES.set(design);
        Ok(DESIGN_BARCODES.get().unwrap())
    }

    pub(super) fn protect_design(index: &mut BarcodeIndex) -> Result<()> {
        let design = design_barcodes()?;
        index.full.retain(|barcode, candidate| {
            matches!(candidate, MatchCandidate::Unique { distance: 0, .. })
                || !design.contains(barcode.as_slice())
        });
        index.shifted.clear();
        Ok(())
    }

    const ANCHOR: &[u8] = b"TTTCTTATATGGGCGTCCGTCGTTGCTCGTAGATGTGTATAAGAGACAG";

    #[derive(Debug, PartialEq, Eq)]
    pub(super) struct Layout<'a> {
        pub barcode: &'a [u8],
        pub umi: &'a [u8],
        pub barcode_start: usize,
        pub trim_start: usize,
    }

    fn mismatch(expected: u8, observed: u8, ct: bool) -> usize {
        usize::from(
            !(expected == observed
                || (ct && expected == b'C' && observed == b'T')
                || (!ct && expected == b'G' && observed == b'A')),
        )
    }

    fn distance(observed: &[u8], ct: bool) -> usize {
        let n = ANCHOR.len();
        if observed.len() == n {
            return ANCHOR
                .iter()
                .zip(observed)
                .map(|(&a, &b)| mismatch(a, b, ct))
                .sum();
        }
        let (short, long, deletion) = if observed.len() < n {
            (observed, ANCHOR, true)
        } else {
            (ANCHOR, observed, false)
        };
        let cost = |i: usize, j: usize| {
            if deletion {
                mismatch(long[j], short[i], ct)
            } else {
                mismatch(short[i], long[j], ct)
            }
        };
        let mut suffix: usize = (0..short.len()).map(|i| cost(i, i + 1)).sum();
        let mut prefix = 0;
        let mut best = 1 + suffix;
        for i in 0..short.len() {
            suffix -= cost(i, i + 1);
            prefix += cost(i, i);
            best = best.min(1 + prefix + suffix);
        }
        best
    }

    pub(super) fn extract(seq: &[u8]) -> std::result::Result<Layout<'_>, TerminalFate> {
        if seq.len() < 28 + ANCHOR.len() - 1 {
            return Err(TerminalFate::Short);
        }
        let mut best = 3;
        let mut boundary = None;
        let mut ambiguous = false;
        for start in 28..=30 {
            for length in 48..=50 {
                let end = start + length;
                if end > seq.len() {
                    continue;
                }
                let score = distance(&seq[start..end], true).min(distance(&seq[start..end], false));
                if score > 2 {
                    continue;
                }
                if score < best {
                    best = score;
                    boundary = Some((start, end));
                    ambiguous = false;
                } else if score == best && boundary != Some((start, end)) {
                    ambiguous = true;
                }
            }
        }
        if ambiguous {
            return Err(TerminalFate::Ambiguous);
        }
        let (start, end) = boundary.ok_or(TerminalFate::Unmatched)?;
        if start < 29 || end + 9 >= seq.len() {
            return Err(TerminalFate::Short);
        }
        let barcode = &seq[start - 29..start - 12];
        let umi = &seq[start - 12..start];
        if !barcode.iter().all(|b| b"ACGT".contains(b)) || !umi.iter().all(|b| b"ACGTN".contains(b))
        {
            return Err(TerminalFate::Unmatched);
        }
        Ok(Layout {
            barcode,
            umi,
            barcode_start: start - 29,
            trim_start: end + 9,
        })
    }

    pub(super) fn match_read<'a>(seq: &'a [u8], index: &'a BarcodeIndex) -> ProcResult<'a> {
        let layout = match extract(seq) {
            Ok(layout) => layout,
            Err(fate) => return ProcResult::rejected(fate, None),
        };
        match index.match_simple(layout.barcode) {
            BarcodeMatch::Unique { id } => {
                ProcResult::assigned(Modality::Dna, id, layout.trim_start, Some(layout.umi))
            }
            BarcodeMatch::Ambiguous => {
                ProcResult::rejected(TerminalFate::Ambiguous, Some(layout.umi))
            }
            BarcodeMatch::None => ProcResult::rejected(TerminalFate::Unmatched, Some(layout.umi)),
        }
    }

    #[derive(Default, Serialize)]
    struct CountMetrics {
        input_reads: u64,
        structured_reads: u64,
        ambiguous: u64,
        unmatched: u64,
        short: u64,
        barcode_with_c: u64,
        umi_with_n: u64,
        unique_barcodes: usize,
    }

    #[derive(Default)]
    struct BarcodeCount {
        count: u64,
        below_q20: [u64; 17],
        q30: [u64; 17],
    }

    pub(super) fn count_barcodes(args: &Args) -> Result<()> {
        anyhow::ensure!(
            args.mode == DemuxMode::DnaOnlyDroplet,
            "--count-only requires dna-only-droplet"
        );
        anyhow::ensure!(args.r2.is_empty(), "--count-only reads R1 only");
        let counts_path = args
            .counts_output
            .as_ref()
            .context("--counts-output is required")?;
        let metrics_path = args
            .count_metrics
            .as_ref()
            .context("--count-metrics is required")?;
        let mut metrics = CountMetrics::default();
        let mut counts: AHashMap<Vec<u8>, BarcodeCount> = AHashMap::new();
        let mut batch = Vec::with_capacity(DEFAULT_CHUNK_SIZE);
        let accumulate = |batch: &mut Vec<(Vec<u8>, Vec<u8>)>,
                          counts: &mut AHashMap<Vec<u8>, BarcodeCount>,
                          metrics: &mut CountMetrics| {
            let layouts: Vec<_> = batch.par_iter().map(|(seq, _)| extract(seq)).collect();
            for (layout, (_, quality)) in layouts.into_iter().zip(batch.iter()) {
                metrics.input_reads += 1;
                match layout {
                    Ok(layout) => {
                        metrics.structured_reads += 1;
                        metrics.barcode_with_c += u64::from(layout.barcode.contains(&b'C'));
                        metrics.umi_with_n += u64::from(layout.umi.contains(&b'N'));
                        let count = counts.entry(layout.barcode.to_vec()).or_default();
                        count.count += 1;
                        for (i, q) in quality[layout.barcode_start..layout.barcode_start + 17]
                            .iter()
                            .enumerate()
                        {
                            count.below_q20[i] += u64::from(*q < 53);
                            count.q30[i] += u64::from(*q >= 63);
                        }
                    }
                    Err(TerminalFate::Ambiguous) => metrics.ambiguous += 1,
                    Err(TerminalFate::Short) => metrics.short += 1,
                    Err(_) => metrics.unmatched += 1,
                }
            }
            batch.clear();
        };
        for path in &args.r1 {
            let mut reader = needletail::parse_fastx_file(path)?;
            while let Some(record) = reader.next() {
                let record = record?;
                let seq = record.seq();
                let qual = record
                    .qual()
                    .context("R1 must contain FASTQ quality scores")?;
                anyhow::ensure!(
                    seq.len() == qual.len(),
                    "FASTQ sequence/quality length mismatch"
                );
                anyhow::ensure!(
                    qual.iter().all(|q| (33..=126).contains(q)),
                    "Invalid FASTQ Phred+33 quality"
                );
                normalize_read_id(record.id())?;
                batch.push((seq.to_vec(), qual.to_vec()));
                if batch.len() == DEFAULT_CHUNK_SIZE {
                    accumulate(&mut batch, &mut counts, &mut metrics);
                }
            }
        }
        accumulate(&mut batch, &mut counts, &mut metrics);
        metrics.unique_barcodes = counts.len();
        let mut rows: Vec<_> = counts.into_iter().collect();
        rows.sort_by(|a, b| b.1.count.cmp(&a.1.count).then(a.0.cmp(&b.0)));
        let mut out = std::io::BufWriter::new(std::fs::File::create(counts_path)?);
        writeln!(out, "barcode\tcount\tbelow_q20\tq30")?;
        for (barcode, count) in rows {
            let list = |values: &[u64; 17]| {
                values
                    .iter()
                    .map(u64::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            };
            writeln!(
                out,
                "{}\t{}\t{}\t{}",
                std::str::from_utf8(&barcode)?,
                count.count,
                list(&count.below_q20),
                list(&count.q30)
            )?;
        }
        out.flush()?;
        std::fs::write(metrics_path, serde_json::to_vec_pretty(&metrics)?)?;
        Ok(())
    }

    const MOLECULE_CAP: usize = 4096;

    fn molecule_key(layout: &Layout<'_>, seq: &[u8]) -> Option<u128> {
        let insert = seq.get(layout.trim_start..layout.trim_start + 32)?;
        layout
            .umi
            .iter()
            .chain(insert)
            .try_fold(0u128, |key, base| {
                let digit = match base {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => return None,
                };
                Some((key << 2) | digit)
            })
    }

    fn molecule_rank(key: u128) -> (u64, u128) {
        let mix = |mut x: u64| {
            x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
            x ^ (x >> 31)
        };
        (mix((key as u64) ^ mix((key >> 64) as u64)), key)
    }

    #[derive(Default)]
    struct Molecules {
        reads: u64,
        saturated: bool,
        keys: std::collections::BTreeSet<(u64, u128)>,
    }

    impl Molecules {
        fn add(&mut self, key: u128) {
            self.reads += 1;
            let ranked = molecule_rank(key);
            if self.keys.len() < MOLECULE_CAP
                || self.keys.last().is_some_and(|last| ranked <= *last)
            {
                self.keys.insert(ranked);
                if self.keys.len() > MOLECULE_CAP {
                    self.keys.pop_last();
                    self.saturated = true;
                }
            } else {
                self.saturated = true;
            }
        }
    }

    fn scan_molecules(
        args: &Args,
        selected: &AHashSet<Vec<u8>>,
        mut consume: impl FnMut(&[u8], u128),
    ) -> Result<()> {
        let mut batch = Vec::with_capacity(DEFAULT_CHUNK_SIZE);
        let mut accumulate = |batch: &mut Vec<Vec<u8>>| {
            let extracted: Vec<_> = batch
                .par_iter()
                .map(|seq| {
                    let layout = extract(seq).ok()?;
                    if !selected.contains(layout.barcode) {
                        return None;
                    }
                    Some((layout.barcode, molecule_key(&layout, seq)?))
                })
                .collect();
            for (barcode, key) in extracted.into_iter().flatten() {
                consume(barcode, key);
            }
            batch.clear();
        };
        for path in &args.r1 {
            let mut reader = needletail::parse_fastx_file(path)?;
            while let Some(record) = reader.next() {
                let record = record?;
                let seq = record.seq();
                anyhow::ensure!(
                    record.qual().is_some_and(|q| q.len() == seq.len()),
                    "Molecule evidence requires FASTQ sequence/quality pairs"
                );
                normalize_read_id(record.id())?;
                batch.push(seq.to_vec());
                if batch.len() == DEFAULT_CHUNK_SIZE {
                    accumulate(&mut batch);
                }
            }
        }
        accumulate(&mut batch);
        Ok(())
    }

    #[derive(serde::Deserialize)]
    #[serde(deny_unknown_fields)]
    struct EvidencePair {
        barcode: String,
        parent: String,
    }

    pub(super) fn molecule_evidence(args: &Args) -> Result<()> {
        anyhow::ensure!(
            args.mode == DemuxMode::DnaOnlyDroplet && args.r2.is_empty() && !args.r1.is_empty(),
            "Droplet molecule evidence reads R1 only"
        );
        let pairs: Vec<EvidencePair> =
            serde_json::from_slice(&std::fs::read(args.evidence_pairs.as_ref().unwrap())?)?;
        let mut children: AHashMap<Vec<u8>, Molecules> = AHashMap::new();
        for pair in &pairs {
            anyhow::ensure!(
                [&pair.barcode, &pair.parent]
                    .into_iter()
                    .all(|bc| bc.len() == 17 && bc.bytes().all(|b| b"AGT".contains(&b)))
                    && pair.barcode != pair.parent,
                "Invalid molecule evidence pair"
            );
            anyhow::ensure!(
                children
                    .insert(pair.barcode.as_bytes().to_vec(), Molecules::default())
                    .is_none(),
                "Duplicate molecule evidence child"
            );
        }
        anyhow::ensure!(
            !pairs.is_empty(),
            "Molecule evidence pairs must not be empty"
        );
        scan_molecules(args, &children.keys().cloned().collect(), |barcode, key| {
            children.get_mut(barcode).unwrap().add(key);
        })?;
        let mut targets: AHashMap<Vec<u8>, AHashMap<u128, Vec<usize>>> = AHashMap::new();
        let mut parent_reads: AHashMap<Vec<u8>, u64> = AHashMap::new();
        let mut shared: Vec<AHashSet<u128>> = pairs.iter().map(|_| AHashSet::new()).collect();
        for (index, pair) in pairs.iter().enumerate() {
            let lookup = targets.entry(pair.parent.as_bytes().to_vec()).or_default();
            parent_reads
                .entry(pair.parent.as_bytes().to_vec())
                .or_default();
            for (_, key) in &children[pair.barcode.as_bytes()].keys {
                lookup.entry(*key).or_default().push(index);
            }
        }
        scan_molecules(args, &targets.keys().cloned().collect(), |barcode, key| {
            *parent_reads.get_mut(barcode).unwrap() += 1;
            if let Some(indices) = targets[barcode].get(&key) {
                for index in indices {
                    shared[*index].insert(key);
                }
            }
        })?;
        let evidence: BTreeMap<_, _> = pairs.iter().enumerate().map(|(index, pair)| {
            let child = &children[pair.barcode.as_bytes()];
            let mut keys: Vec<_> = shared[index].iter().copied().collect();
            keys.sort_unstable();
            (pair.barcode.clone(), serde_json::json!({
                "parent": pair.parent,
                "child_eligible_reads": child.reads,
                "child_sampled_signatures": child.keys.len(),
                "child_saturated": child.saturated,
                "parent_eligible_reads": parent_reads[pair.parent.as_bytes()],
                "shared_signatures": keys.into_iter().map(|key| format!("{key:022x}")).collect::<Vec<_>>()
            }))
        }).collect();
        std::fs::write(
            args.evidence_output.as_ref().unwrap(),
            serde_json::to_vec(
                &serde_json::json!({"signature": "exact_umi12_insert32", "child_cap": MOLECULE_CAP,
                    "parent_scan": "complete", "pairs": evidence}),
            )?,
        )?;
        Ok(())
    }

    #[test]
    fn molecule_sample_is_bounded_order_independent_and_duplicate_safe() {
        let mut forward = Molecules::default();
        let mut reverse = Molecules::default();
        for key in 0..10000u128 {
            forward.add(key);
            forward.add(key);
        }
        for key in (0..10000u128).rev() {
            reverse.add(key);
        }
        assert_eq!(forward.keys, reverse.keys);
        assert_eq!(forward.keys.len(), MOLECULE_CAP);
        assert!(forward.saturated && reverse.saturated);
        assert_eq!(forward.reads, 20000);
        let mut small = Molecules::default();
        for _ in 0..10000 {
            small.add(123);
        }
        assert_eq!(small.keys.len(), 1);
        assert!(!small.saturated);
    }

    #[test]
    fn molecule_signature_uses_complete_umi_and_insert_without_ambiguous_bases() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut layout = Layout {
            barcode: b"AAAAAAAAAAAAAAAAA",
            umi: b"ACGTACGTACGT",
            barcode_start: 0,
            trim_start: 0,
        };
        let key = molecule_key(&layout, seq).unwrap();
        layout.umi = b"TCGTACGTACGT";
        assert_ne!(molecule_key(&layout, seq).unwrap(), key);
        layout.umi = b"NCGTACGTACGT";
        assert!(molecule_key(&layout, seq).is_none());
        layout.umi = b"ACGTACGTACGT";
        assert!(molecule_key(&layout, &seq[..31]).is_none());
        assert!(molecule_key(&layout, b"NCGTACGTACGTACGTACGTACGTACGTACGT").is_none());
    }
}

use barcode_map::*;
use fastq_io::*;
use matcher::*;
use output::*;

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use clap::Parser;
use cli::{Args, DemuxMode};
use crossbeam_channel::{bounded, Receiver, Sender};
use libdeflate::{CompressionLvl, Compressor};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::BTreeMap;
use std::hash::Hash;
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};
use std::thread::JoinHandle;

static DNA_NATIVE_TABLE: OnceLock<[u8; 256]> = OnceLock::new();
static DNA_BISULFITE_TABLE: OnceLock<[u8; 256]> = OnceLock::new();
static RNA_BISULFITE_TABLE: OnceLock<[u8; 256]> = OnceLock::new();

const ME_SEQ: &[u8] = b"AGATGTGTATAAGAGACAG";
const TSO_SEQ: &[u8] = b"GTCTAACGCGTTAC";
const LEN_UMI: usize = 8;
const LEN_SPACER: usize = 5;
const LEN_GAP_AFTER_ME: usize = 9;
const MIN_R1_LEN: usize = 20;
const DEFAULT_CHUNK_SIZE: usize = 20_000;
const PROGRESS_EVERY_READS: u64 = 5_000_000;
const BUILD_SOURCE_REVISION: &str = match option_env!("DNA_PIPELINE_SOURCE_REVISION") {
    Some(value) => value,
    None => "unrecorded",
};

#[derive(Debug, Clone)]
struct Config {
    r1: Vec<PathBuf>,
    r2: Vec<PathBuf>,
    barcode_map: PathBuf,
    dna_out_dir: PathBuf,
    rna_out_dir: PathBuf,
    json_report: PathBuf,
    sample_name: String,
    min_matched_read_pairs: u64,
    threads: usize,
    dna_w_spacer_len: usize,
    mode: DemuxMode,
    chunk_size: usize,
}

impl Config {
    fn from_args(args: Args) -> Result<Self> {
        if args.sample_name.is_empty()
            || matches!(args.sample_name.as_str(), "." | "..")
            || !args
                .sample_name
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'.' | b'-'))
        {
            anyhow::bail!("--sample-name must be a safe non-empty path component");
        }
        if args.threads == 0 {
            anyhow::bail!("--threads must be greater than 0");
        }
        if args.dna_w_spacer_len > 20 {
            anyhow::bail!(
                "--dna-w-spacer-len is unexpectedly large: {}",
                args.dna_w_spacer_len
            );
        }
        if args.r1.is_empty() || args.r2.is_empty() {
            anyhow::bail!("At least one --r1/--r2 FASTQ pair is required");
        }
        if args.r1.len() != args.r2.len() {
            anyhow::bail!(
                "Multi-lane input mismatch: received {} --r1 file(s) but {} --r2 file(s)",
                args.r1.len(),
                args.r2.len()
            );
        }
        for path in &args.r1 {
            ensure_gzip_fastq(path, "--r1")?;
        }
        for path in &args.r2 {
            ensure_gzip_fastq(path, "--r2")?;
        }
        let barcode_map = args
            .barcode_map
            .context("--barcode-map is required for demux")?;
        let dna_out_dir = args
            .dna_out_dir
            .context("--dna-out-dir is required for demux")?;
        let rna_out_dir = args
            .rna_out_dir
            .context("--rna-out-dir is required for demux")?;
        let json_report = args
            .json_report
            .context("--json-report is required for demux")?;
        ensure_readable_file(&barcode_map, "--barcode-map")?;
        ensure_parent_dir(&json_report, "--json-report")?;
        Ok(Self {
            r1: args.r1,
            r2: args.r2,
            barcode_map,
            dna_out_dir,
            rna_out_dir,
            json_report,
            sample_name: args.sample_name,
            min_matched_read_pairs: args.min_matched_read_pairs,
            threads: args.threads,
            dna_w_spacer_len: args.dna_w_spacer_len,
            mode: args.mode,
            chunk_size: DEFAULT_CHUNK_SIZE,
        })
    }
}

#[derive(Debug, Clone)]
struct CellRow {
    dna_barcode: String,
    rna_barcode: String,
    plate_id: String,
    cell_order: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Modality {
    Dna,
    Rna,
}

impl Modality {
    fn as_str(self) -> &'static str {
        match self {
            Self::Dna => "DNA",
            Self::Rna => "RNA",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum BucketKey {
    Dna(String),
    Rna(String),
    Unmatched,
}

impl BucketKey {
    fn is_matched(&self) -> bool {
        matches!(self, Self::Dna(_) | Self::Rna(_))
    }

    fn label(&self) -> &str {
        match self {
            Self::Dna(bc) | Self::Rna(bc) => bc.as_str(),
            Self::Unmatched => "unmatched",
        }
    }

    fn cell_sample_id<'a>(&self, paths: &'a OutputPaths) -> Result<&'a str> {
        match self {
            Self::Dna(barcode) => paths
                .dna_cell_ids
                .get(barcode)
                .map(String::as_str)
                .with_context(|| format!("No canonical Sample ID for DNA barcode {}", barcode)),
            Self::Rna(barcode) => paths
                .rna_cell_ids
                .get(barcode)
                .map(String::as_str)
                .with_context(|| format!("No canonical Sample ID for RNA barcode {}", barcode)),
            Self::Unmatched => Ok(paths.sample_name.as_str()),
        }
    }

    fn output_dir(&self, paths: &OutputPaths) -> Result<PathBuf> {
        match self {
            Self::Dna(_) => Ok(paths.dna_out_dir.join(self.cell_sample_id(paths)?)),
            Self::Rna(_) => Ok(paths.rna_out_dir.join(self.cell_sample_id(paths)?)),
            Self::Unmatched => Ok(paths
                .dna_out_dir
                .join(format!("{}_unmatched", paths.sample_name))),
        }
    }

    fn file_stem(&self, paths: &OutputPaths) -> Result<String> {
        match self {
            Self::Dna(_) | Self::Rna(_) => Ok(self.cell_sample_id(paths)?.to_string()),
            Self::Unmatched => Ok(format!("{}_unmatched", paths.sample_name)),
        }
    }
}

#[derive(Debug, Clone)]
struct OutputPaths {
    dna_out_dir: PathBuf,
    rna_out_dir: PathBuf,
    json_report: PathBuf,
    sample_name: String,
    dna_cell_ids: AHashMap<String, String>,
    rna_cell_ids: AHashMap<String, String>,
    min_matched_read_pairs: u64,
    mode: DemuxMode,
    input_fastq_pairs: usize,
}

#[derive(Debug, Default, Clone, Copy, Serialize)]
struct FateCounts {
    dna_assigned: u64,
    rna_assigned: u64,
    ambiguous: u64,
    unmatched: u64,
    short: u64,
}

impl FateCounts {
    fn record(&mut self, fate: TerminalFate) {
        match fate {
            TerminalFate::DnaAssigned => self.dna_assigned += 1,
            TerminalFate::RnaAssigned => self.rna_assigned += 1,
            TerminalFate::Ambiguous => self.ambiguous += 1,
            TerminalFate::Unmatched => self.unmatched += 1,
            TerminalFate::Short => self.short += 1,
        }
    }

    fn add_assign(&mut self, other: Self) {
        self.dna_assigned += other.dna_assigned;
        self.rna_assigned += other.rna_assigned;
        self.ambiguous += other.ambiguous;
        self.unmatched += other.unmatched;
        self.short += other.short;
    }

    fn input_read_pairs(self) -> u64 {
        self.dna_assigned + self.rna_assigned + self.ambiguous + self.unmatched + self.short
    }

    fn usable_read_pairs(self) -> u64 {
        self.input_read_pairs().saturating_sub(self.short)
    }

    fn matched_read_pairs(self) -> u64 {
        self.dna_assigned + self.rna_assigned
    }
}

#[derive(Serialize)]
struct DemuxReport {
    schema_version: u32,
    build_version: String,
    source_revision: String,
    mode: String,
    retention_policy: String,
    retention_threshold: u64,
    input_fastq_pairs: usize,
    input_read_pairs: u64,
    usable_read_pairs: u64,
    matched_reads: u64,
    read_fates: FateCounts,
    samples: Vec<SampleStats>,
}

#[derive(Serialize)]
struct SampleStats {
    sample_name: String,
    cell_sample_id: String,
    plate_id: String,
    cell_order: String,
    dna_barcode: String,
    rna_barcode: String,
    dna_read_count: u64,
    rna_read_count: u64,
    dna_status: String,
    rna_status: String,
}

struct InputChunk {
    r1_seqs: Vec<u8>,
    r1_quals: Vec<u8>,
    r1_ids: Vec<u8>,
    r1_offsets: Vec<usize>,
    r1_id_offsets: Vec<usize>,
    r2_seqs: Vec<u8>,
    r2_quals: Vec<u8>,
    r2_ids: Vec<u8>,
    r2_offsets: Vec<usize>,
    r2_id_offsets: Vec<usize>,
    count: usize,
}

impl InputChunk {
    fn new(capacity: usize) -> Self {
        Self {
            r1_seqs: Vec::with_capacity(capacity * 150),
            r1_quals: Vec::with_capacity(capacity * 150),
            r1_ids: Vec::with_capacity(capacity * 50),
            r1_offsets: vec![0],
            r1_id_offsets: vec![0],
            r2_seqs: Vec::with_capacity(capacity * 150),
            r2_quals: Vec::with_capacity(capacity * 150),
            r2_ids: Vec::with_capacity(capacity * 50),
            r2_offsets: vec![0],
            r2_id_offsets: vec![0],
            count: 0,
        }
    }

    fn clear(&mut self) {
        self.r1_seqs.clear();
        self.r1_quals.clear();
        self.r1_ids.clear();
        self.r1_offsets.truncate(1);
        self.r1_id_offsets.truncate(1);
        self.r2_seqs.clear();
        self.r2_quals.clear();
        self.r2_ids.clear();
        self.r2_offsets.truncate(1);
        self.r2_id_offsets.truncate(1);
        self.count = 0;
    }

    fn push_pair(
        &mut self,
        r1_id: &[u8],
        r1_seq: &[u8],
        r1_qual: &[u8],
        r2_id: &[u8],
        r2_seq: &[u8],
        r2_qual: &[u8],
    ) {
        self.r1_ids.extend_from_slice(r1_id);
        self.r1_seqs.extend_from_slice(r1_seq);
        self.r1_quals.extend_from_slice(r1_qual);
        self.r1_id_offsets.push(self.r1_ids.len());
        self.r1_offsets.push(self.r1_seqs.len());

        self.r2_ids.extend_from_slice(r2_id);
        self.r2_seqs.extend_from_slice(r2_seq);
        self.r2_quals.extend_from_slice(r2_qual);
        self.r2_id_offsets.push(self.r2_ids.len());
        self.r2_offsets.push(self.r2_seqs.len());

        self.count += 1;
    }
}

struct FastqBucket {
    r1_raw: Vec<u8>,
    r2_raw: Vec<u8>,
    r1_gzip: Vec<u8>,
    r2_gzip: Vec<u8>,
    count: u64,
}

impl FastqBucket {
    fn new() -> Self {
        Self {
            r1_raw: Vec::new(),
            r2_raw: Vec::new(),
            r1_gzip: Vec::new(),
            r2_gzip: Vec::new(),
            count: 0,
        }
    }

    fn clear(&mut self) {
        self.r1_raw.clear();
        self.r2_raw.clear();
        self.r1_gzip.clear();
        self.r2_gzip.clear();
        self.count = 0;
    }
}

struct AccumulationBuffer {
    map: AHashMap<BucketKey, FastqBucket>,
    fates: FateCounts,
}

type WorkChunk = (u64, InputChunk, AccumulationBuffer);
type CompletedChunk = (u64, AccumulationBuffer);

fn insert_completed_chunk<T>(
    pending: &mut BTreeMap<u64, T>,
    next_chunk_id: &mut u64,
    chunk_id: u64,
    item: T,
) -> Result<Vec<T>> {
    if chunk_id < *next_chunk_id || pending.contains_key(&chunk_id) {
        anyhow::bail!("Worker produced duplicate or stale chunk id {chunk_id}");
    }
    pending.insert(chunk_id, item);

    let mut ready = Vec::new();
    while let Some(item) = pending.remove(&*next_chunk_id) {
        ready.push(item);
        *next_chunk_id = (*next_chunk_id)
            .checked_add(1)
            .context("Completed chunk id overflow")?;
    }
    Ok(ready)
}

impl AccumulationBuffer {
    fn new() -> Self {
        Self {
            map: AHashMap::with_capacity(128),
            fates: FateCounts::default(),
        }
    }

    fn bucket_mut(&mut self, key: BucketKey) -> &mut FastqBucket {
        self.map.entry(key).or_insert_with(FastqBucket::new)
    }

    fn clear(&mut self) {
        for bucket in self.map.values_mut() {
            bucket.clear();
        }
        self.fates = FateCounts::default();
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MatchCandidate {
    Unique { id: String, distance: u8 },
    Ambiguous { distance: u8 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BarcodeMatch<'a> {
    Unique { id: &'a str },
    Ambiguous,
    None,
}

type DistMap = AHashMap<Vec<u8>, MatchCandidate>;

struct BarcodeIndex {
    full: DistMap,
    shifted: DistMap,
    barcode_len: usize,
    modality: Modality,
}

impl BarcodeIndex {
    fn empty(modality: Modality) -> Self {
        Self {
            full: DistMap::new(),
            shifted: DistMap::new(),
            barcode_len: 0,
            modality,
        }
    }

    fn from_barcodes(
        barcodes: impl IntoIterator<Item = String>,
        barcode_len: usize,
        modality: Modality,
    ) -> Self {
        let mut full = DistMap::new();
        let mut shifted = DistMap::new();
        let bases = b"ACGT";

        for barcode in barcodes {
            let seq = barcode.trim().to_ascii_uppercase().into_bytes();
            if seq.len() < barcode_len {
                continue;
            }
            let id = String::from_utf8_lossy(&seq).to_string();
            let key = seq[..barcode_len].to_vec();

            Self::insert_candidate(&mut full, key.clone(), &id, 0);
            Self::gen_one_substitution_candidates(&key, bases, |candidate, dist| {
                Self::insert_candidate(&mut full, candidate, &id, dist);
            });

            if modality == Modality::Dna && barcode_len > 1 {
                let shifted_key = seq[1..barcode_len].to_vec();
                Self::insert_candidate(&mut shifted, shifted_key, &id, 0);
            }
        }

        Self {
            full,
            shifted,
            barcode_len,
            modality,
        }
    }

    fn len(&self) -> usize {
        self.barcode_len
    }

    fn match_cascade(&self, seq: &[u8]) -> BarcodeMatch<'_> {
        if seq
            .iter()
            .any(|base| !matches!(*base, b'A' | b'C' | b'G' | b'T'))
        {
            return BarcodeMatch::None;
        }
        let direct = Self::match_candidate(self.full.get(seq));
        if direct != BarcodeMatch::None {
            return direct;
        }
        if self.modality == Modality::Dna && seq.len() + 1 == self.barcode_len {
            return Self::match_candidate(self.shifted.get(seq));
        }
        BarcodeMatch::None
    }

    fn match_simple(&self, seq: &[u8]) -> BarcodeMatch<'_> {
        if seq
            .iter()
            .any(|base| !matches!(*base, b'A' | b'C' | b'G' | b'T'))
        {
            return BarcodeMatch::None;
        }
        Self::match_candidate(self.full.get(seq))
    }

    fn match_candidate(candidate: Option<&MatchCandidate>) -> BarcodeMatch<'_> {
        match candidate {
            Some(MatchCandidate::Unique { id, .. }) => BarcodeMatch::Unique { id: id.as_str() },
            Some(MatchCandidate::Ambiguous { .. }) => BarcodeMatch::Ambiguous,
            None => BarcodeMatch::None,
        }
    }

    fn insert_candidate(map: &mut DistMap, seq: Vec<u8>, id: &str, dist: u8) {
        use std::collections::hash_map::Entry;

        match map.entry(seq) {
            Entry::Vacant(entry) => {
                entry.insert(MatchCandidate::Unique {
                    id: id.to_string(),
                    distance: dist,
                });
            }
            Entry::Occupied(mut entry) => {
                let replacement = match entry.get() {
                    MatchCandidate::Unique {
                        id: existing_id,
                        distance: existing_dist,
                    } if dist < *existing_dist => Some(MatchCandidate::Unique {
                        id: id.to_string(),
                        distance: dist,
                    }),
                    MatchCandidate::Unique {
                        id: existing_id,
                        distance: existing_dist,
                    } if dist == *existing_dist && existing_id != id => {
                        Some(MatchCandidate::Ambiguous { distance: dist })
                    }
                    MatchCandidate::Ambiguous {
                        distance: existing_dist,
                    } if dist < *existing_dist => Some(MatchCandidate::Unique {
                        id: id.to_string(),
                        distance: dist,
                    }),
                    _ => None,
                };
                if let Some(candidate) = replacement {
                    entry.insert(candidate);
                }
            }
        }
    }

    fn gen_one_substitution_candidates<F>(seq: &[u8], bases: &[u8], mut insert: F)
    where
        F: FnMut(Vec<u8>, u8),
    {
        for pos in 0..seq.len() {
            for &base in bases {
                if base == seq[pos] {
                    continue;
                }
                let mut candidate = seq.to_vec();
                candidate[pos] = base;
                insert(candidate, 1);
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TerminalFate {
    DnaAssigned,
    RnaAssigned,
    Ambiguous,
    Unmatched,
    Short,
}

struct ProcResult<'a> {
    trim_start: usize,
    category: &'a str,
    modality: Option<Modality>,
    umi: Option<&'a [u8]>,
    is_valid: bool,
    fate: TerminalFate,
}

impl<'a> ProcResult<'a> {
    fn unmatched() -> Self {
        Self {
            trim_start: 0,
            category: "unmatched",
            modality: None,
            umi: None,
            is_valid: false,
            fate: TerminalFate::Unmatched,
        }
    }

    fn cross_layout_ambiguous() -> Self {
        Self {
            fate: TerminalFate::Ambiguous,
            ..Self::unmatched()
        }
    }

    fn rejected(fate: TerminalFate, umi: Option<&'a [u8]>) -> Self {
        Self {
            fate,
            umi,
            ..Self::unmatched()
        }
    }

    fn assigned(
        modality: Modality,
        assigned_barcode: &'a str,
        trim_start: usize,
        umi: Option<&'a [u8]>,
    ) -> Self {
        Self {
            trim_start,
            category: assigned_barcode,
            modality: Some(modality),
            umi,
            is_valid: true,
            fate: match modality {
                Modality::Dna => TerminalFate::DnaAssigned,
                Modality::Rna => TerminalFate::RnaAssigned,
            },
        }
    }

    fn output_umi(&self) -> Option<&'a [u8]> {
        if self.is_valid {
            self.umi
        } else {
            None
        }
    }
}

fn main() {
    if std::env::args_os().any(|argument| argument == "--build-info") {
        println!(
            "{}",
            serde_json::json!({
                "name": env!("CARGO_PKG_NAME"),
                "version": env!("CARGO_PKG_VERSION"),
                "source_revision": BUILD_SOURCE_REVISION,
            })
        );
        return;
    }
    env_logger::init();
    if let Err(err) = run() {
        eprintln!("{err:#}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    init_tables();

    let args = Args::parse();
    if args.threads == 0 {
        anyhow::bail!("--threads must be greater than 0");
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .context("Cannot initialize Rayon global thread pool")?;
    if args.count_only {
        return droplet::count_barcodes(&args);
    }
    if args.evidence_pairs.is_some() {
        return droplet::molecule_evidence(&args);
    }
    let config = Config::from_args(args)?;

    log::info!(
        "Start demux. Sample: {}; build_version={}; source_revision={}",
        config.sample_name,
        env!("CARGO_PKG_VERSION"),
        BUILD_SOURCE_REVISION
    );

    let cell_rows = load_barcode_map(&config.barcode_map, config.mode)?;
    let (dna_bc_lengths, rna_bc_len) = resolve_barcode_lengths(&cell_rows, config.mode)?;
    validate_effective_barcode_prefixes(&cell_rows, &dna_bc_lengths, rna_bc_len, config.mode)?;
    prepare_output_roots(&config)?;
    let dna_bc_lengths_label = dna_bc_lengths
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(",");
    if config.mode == DemuxMode::DnaRna {
        log::info!(
            "Using barcode lengths: DNA={}bp, RNA={}bp; DNA barcode-to-ME spacer={}bp",
            dna_bc_lengths_label,
            rna_bc_len,
            config.dna_w_spacer_len
        );
    } else {
        log::info!(
            "Using barcode lengths: DNA={}bp; RNA demux disabled; DNA barcode-to-ME spacer={}bp",
            dna_bc_lengths_label,
            config.dna_w_spacer_len
        );
    }
    log::info!("Using demultiplexing mode: {:?}", config.mode);

    let mixed_dna_lengths = dna_bc_lengths.len() > 1;
    let dna_indices = Arc::new(
        dna_bc_lengths
            .iter()
            .map(|&dna_bc_len| {
                let dna_barcodes = cell_rows
                    .iter()
                    .filter(move |row| !mixed_dna_lengths || row.dna_barcode.len() == dna_bc_len)
                    .map(|row| row.dna_barcode.clone());
                let mut index =
                    BarcodeIndex::from_barcodes(dna_barcodes, dna_bc_len, Modality::Dna);
                if config.mode == DemuxMode::DnaOnlyDroplet {
                    droplet::protect_design(&mut index)?;
                }
                Ok(index)
            })
            .collect::<Result<Vec<_>>>()?,
    );
    let rna_idx = if config.mode == DemuxMode::DnaRna {
        let rna_barcodes = cell_rows.iter().map(|row| row.rna_barcode.clone());
        Arc::new(BarcodeIndex::from_barcodes(
            rna_barcodes,
            rna_bc_len,
            Modality::Rna,
        ))
    } else {
        Arc::new(BarcodeIndex::empty(Modality::Rna))
    };

    let dna_cell_ids = cell_rows
        .iter()
        .map(|row| {
            (
                row.dna_barcode.clone(),
                format!(
                    "{}_{}",
                    config.sample_name,
                    if config.mode == DemuxMode::DnaOnlyDroplet {
                        &row.dna_barcode
                    } else {
                        &row.plate_id
                    }
                ),
            )
        })
        .collect::<AHashMap<_, _>>();
    let rna_cell_ids = cell_rows
        .iter()
        .filter(|row| !row.rna_barcode.is_empty())
        .map(|row| {
            (
                row.rna_barcode.clone(),
                format!("{}_{}", config.sample_name, row.plate_id),
            )
        })
        .collect::<AHashMap<_, _>>();

    let paths = OutputPaths {
        dna_out_dir: config.dna_out_dir.clone(),
        rna_out_dir: config.rna_out_dir.clone(),
        json_report: config.json_report.clone(),
        sample_name: config.sample_name.clone(),
        dna_cell_ids,
        rna_cell_ids,
        min_matched_read_pairs: config.min_matched_read_pairs,
        mode: config.mode,
        input_fastq_pairs: config.r1.len(),
    };

    run_pipeline(config, cell_rows, dna_indices, rna_idx, paths)
}

fn run_pipeline(
    config: Config,
    cell_rows: Vec<CellRow>,
    dna_indices: Arc<Vec<BarcodeIndex>>,
    rna_idx: Arc<BarcodeIndex>,
    paths: OutputPaths,
) -> Result<()> {
    let pool_capacity = config
        .threads
        .checked_mul(2)
        .context("--threads is too large: channel capacity overflow")?
        .max(1);

    let (tx_input_recycle, rx_input_recycle) = bounded::<InputChunk>(pool_capacity);
    let (tx_acc_recycle, rx_acc_recycle) = bounded::<AccumulationBuffer>(pool_capacity);
    let (tx_work, rx_work) = bounded::<WorkChunk>(pool_capacity);
    let (tx_write, rx_write) = bounded::<CompletedChunk>(pool_capacity);

    for _ in 0..pool_capacity {
        tx_input_recycle
            .send(InputChunk::new(config.chunk_size))
            .context("Cannot initialize input chunk pool")?;
        tx_acc_recycle
            .send(AccumulationBuffer::new())
            .context("Cannot initialize accumulation buffer pool")?;
    }

    let dna_lut = if config.mode == DemuxMode::DnaOnlyTaps {
        DNA_NATIVE_TABLE.get()
    } else {
        DNA_BISULFITE_TABLE.get()
    }
    .context("DNA table was not initialized")?;
    let rna_lut = RNA_BISULFITE_TABLE
        .get()
        .context("RNA bisulfite table was not initialized")?;

    let worker_handle = spawn_worker(
        rx_work,
        tx_write,
        tx_input_recycle,
        dna_indices,
        rna_idx,
        dna_lut,
        rna_lut,
        config.dna_w_spacer_len,
        config.mode,
    );

    let writer_handle = spawn_writer(rx_write, tx_acc_recycle, paths, cell_rows);

    let read_result = read_fastq_pairs(
        &config.r1,
        &config.r2,
        config.chunk_size,
        &rx_input_recycle,
        &rx_acc_recycle,
        &tx_work,
    );
    drop(tx_work);

    if read_result.is_err() {
        drop(rx_input_recycle);
        drop(rx_acc_recycle);
    }

    let worker_result = match join_thread(worker_handle, "worker") {
        Ok(result) => result,
        Err(err) => Err(err),
    };
    let writer_result = match join_thread(writer_handle, "writer") {
        Ok(result) => result,
        Err(err) => Err(err),
    };

    if let Err(writer_err) = writer_result.as_ref() {
        let mut message = format!("Demux writer failed: {writer_err:#}");
        if let Err(worker_err) = worker_result.as_ref() {
            message.push_str(&format!("\nDemux worker also reported: {worker_err:#}"));
        }
        if let Err(reader_err) = read_result.as_ref() {
            message.push_str(&format!("\nFASTQ reader also reported: {reader_err:#}"));
        }
        anyhow::bail!(message);
    }
    if let Err(worker_err) = worker_result.as_ref() {
        let mut message = format!("Demux worker failed: {worker_err:#}");
        if let Err(reader_err) = read_result.as_ref() {
            message.push_str(&format!("\nFASTQ reader also reported: {reader_err:#}"));
        }
        anyhow::bail!(message);
    }
    read_result?;

    let matched_reads = writer_result.expect("writer_result was checked above");
    log::info!("Done. Matched reads: {}", matched_reads);
    Ok(())
}

fn spawn_worker(
    rx_work: Receiver<WorkChunk>,
    tx_write: Sender<CompletedChunk>,
    tx_input_recycle: Sender<InputChunk>,
    dna_indices: Arc<Vec<BarcodeIndex>>,
    rna_idx: Arc<BarcodeIndex>,
    dna_lut: &'static [u8; 256],
    rna_lut: &'static [u8; 256],
    dna_w_spacer_len: usize,
    demux_mode: DemuxMode,
) -> JoinHandle<Result<()>> {
    std::thread::spawn(move || {
        rx_work.into_iter().par_bridge().try_for_each_with(
            (
                tx_write,
                tx_input_recycle,
                dna_indices,
                rna_idx,
                dna_w_spacer_len,
                demux_mode,
            ),
            |(writer, recycle_input, dna_indexes, rna_index, spacer_len, mode),
             (chunk_id, mut input_chunk, mut acc_buf)| {
                process_chunk(
                    &input_chunk,
                    &mut acc_buf,
                    dna_indexes,
                    rna_index,
                    dna_lut,
                    rna_lut,
                    *spacer_len,
                    *mode,
                )?;
                compress_accumulation(&mut acc_buf)?;
                input_chunk.clear();
                recycle_input
                    .send(input_chunk)
                    .context("Input chunk recycle channel is closed")?;
                writer
                    .send((chunk_id, acc_buf))
                    .context("Writer channel is closed")?;
                Ok::<(), anyhow::Error>(())
            },
        )?;
        Ok(())
    })
}

fn process_chunk(
    input: &InputChunk,
    acc: &mut AccumulationBuffer,
    dna_indices: &[BarcodeIndex],
    rna_idx: &BarcodeIndex,
    dna_lut: &'static [u8; 256],
    rna_lut: &'static [u8; 256],
    dna_w_spacer_len: usize,
    demux_mode: DemuxMode,
) -> Result<()> {
    for i in 0..input.count {
        let r1_start = input.r1_offsets[i];
        let r1_end = input.r1_offsets[i + 1];
        let r2_start = input.r2_offsets[i];
        let r2_end = input.r2_offsets[i + 1];
        let r1_id_start = input.r1_id_offsets[i];
        let r1_id_end = input.r1_id_offsets[i + 1];
        let r2_id_start = input.r2_id_offsets[i];
        let r2_id_end = input.r2_id_offsets[i + 1];

        let r1_seq = &input.r1_seqs[r1_start..r1_end];
        if r1_seq.len() < MIN_R1_LEN {
            acc.fates.record(TerminalFate::Short);
            continue;
        }
        let mut result = process_read(
            r1_seq,
            dna_indices,
            rna_idx,
            dna_lut,
            rna_lut,
            dna_w_spacer_len,
            demux_mode,
        );
        if demux_mode == DemuxMode::DnaOnlyDroplet && result.is_valid && r2_end - r2_start <= 9 {
            result = ProcResult::rejected(TerminalFate::Short, result.umi);
        }
        acc.fates.record(result.fate);
        if result.fate == TerminalFate::Short {
            continue;
        }
        let key = bucket_key_for_result(&result);
        let is_matched_bucket = key.is_matched();
        let bucket = acc.bucket_mut(key);
        bucket.count += 1;

        if !is_matched_bucket {
            continue;
        }

        write_fastq_record(
            &mut bucket.r1_raw,
            if demux_mode == DemuxMode::DnaOnlyDroplet {
                normalize_read_id(&input.r1_ids[r1_id_start..r1_id_end])?
            } else {
                &input.r1_ids[r1_id_start..r1_id_end]
            },
            r1_seq,
            &input.r1_quals[r1_start..r1_end],
            result.trim_start,
            result.output_umi(),
            result.is_valid,
        )?;
        write_fastq_record(
            &mut bucket.r2_raw,
            if demux_mode == DemuxMode::DnaOnlyDroplet {
                normalize_read_id(&input.r2_ids[r2_id_start..r2_id_end])?
            } else {
                &input.r2_ids[r2_id_start..r2_id_end]
            },
            &input.r2_seqs[r2_start..r2_end],
            &input.r2_quals[r2_start..r2_end],
            if demux_mode == DemuxMode::DnaOnlyDroplet {
                9
            } else {
                0
            },
            result.output_umi(),
            demux_mode == DemuxMode::DnaOnlyDroplet,
        )?;
    }
    Ok(())
}

fn bucket_key_for_result(result: &ProcResult<'_>) -> BucketKey {
    if result.is_valid {
        match result.modality {
            Some(Modality::Dna) => BucketKey::Dna(result.category.to_string()),
            Some(Modality::Rna) => BucketKey::Rna(result.category.to_string()),
            None => BucketKey::Unmatched,
        }
    } else {
        BucketKey::Unmatched
    }
}

fn init_tables() {
    DNA_NATIVE_TABLE.get_or_init(|| {
        let mut table = [0u8; 256];
        for (idx, value) in table.iter_mut().enumerate() {
            *value = (idx as u8).to_ascii_uppercase();
        }
        table
    });
    DNA_BISULFITE_TABLE.get_or_init(|| {
        let mut table = [0u8; 256];
        for (idx, value) in table.iter_mut().enumerate() {
            *value = idx as u8;
        }
        table[b'G' as usize] = b'A';
        table[b'g' as usize] = b'A';
        table
    });
    RNA_BISULFITE_TABLE.get_or_init(|| {
        let mut table = [0u8; 256];
        for (idx, value) in table.iter_mut().enumerate() {
            *value = idx as u8;
        }
        table[b'C' as usize] = b'T';
        table[b'c' as usize] = b'T';
        table
    });
}

fn ensure_gzip_fastq(path: &Path, label: &str) -> Result<()> {
    ensure_readable_file(path, label)?;
    let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_ascii_lowercase();
    if !name.ends_with(".fastq.gz") && !name.ends_with(".fq.gz") {
        anyhow::bail!(
            "{} requires a .fastq.gz or .fq.gz suffix: {}",
            label,
            path.display()
        );
    }
    let mut file = std::fs::File::open(path)
        .with_context(|| format!("Cannot open {} gzip FASTQ: {}", label, path.display()))?;
    let mut magic = [0u8; 2];
    std::io::Read::read_exact(&mut file, &mut magic)
        .with_context(|| format!("Cannot read {} gzip header: {}", label, path.display()))?;
    if magic != [0x1f, 0x8b] {
        anyhow::bail!(
            "{} is not gzip FASTQ (expected gzip magic 0x1f 0x8b): {}",
            label,
            path.display()
        );
    }
    Ok(())
}

fn ensure_readable_file(path: &Path, label: &str) -> Result<()> {
    let metadata = std::fs::metadata(path).with_context(|| {
        format!(
            "{} does not exist or is not accessible: {}",
            label,
            path.display()
        )
    })?;
    if !metadata.is_file() {
        anyhow::bail!("{} is not a file: {}", label, path.display());
    }
    Ok(())
}

fn join_thread<T>(handle: JoinHandle<T>, name: &str) -> Result<T> {
    handle.join().map_err(|panic_payload| {
        let message = if let Some(text) = panic_payload.downcast_ref::<&str>() {
            (*text).to_string()
        } else if let Some(text) = panic_payload.downcast_ref::<String>() {
            text.clone()
        } else {
            "unknown panic payload".to_string()
        };
        anyhow::anyhow!("{} thread panicked: {}", name, message)
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn droplet_profiles_offsets_and_edits() {
        let anchor = b"TTTCTTATATGGGCGTCCGTCGTTGCTCGTAGATGTGTATAAGAGACAG";
        let cb = b"AGTAGTAGTAGTAGTAG";
        let umi = b"ACGTNACGTACG";
        for ct in [true, false] {
            let converted: Vec<u8> = anchor
                .iter()
                .map(|&b| {
                    if ct && b == b'C' {
                        b'T'
                    } else if !ct && b == b'G' {
                        b'A'
                    } else {
                        b
                    }
                })
                .collect();
            let mut seq = [
                cb.as_slice(),
                umi,
                &converted,
                b"AAAAAAAAA",
                b"GCGCGCGCGCGCGCGCGCGC",
            ]
            .concat();
            let found = droplet::extract(&seq).unwrap();
            assert_eq!(found.barcode, cb);
            assert_eq!(found.umi, umi);
            assert_eq!(found.trim_start, 87);
            seq.insert(0, b'G');
            assert_eq!(droplet::extract(&seq).unwrap().trim_start, 88);
            seq.remove(0);
            seq.remove(0);
            assert!(droplet::extract(&seq).is_err());
        }
        let original = [
            cb.as_slice(),
            umi,
            anchor,
            b"AAAAAAAAA",
            b"GCGCGCGCGCGCGCGCGCGC",
        ]
        .concat();
        let mut inserted = original.clone();
        inserted.insert(51, b'A');
        assert_eq!(droplet::extract(&inserted).unwrap().trim_start, 88);
        let mut deleted = original.clone();
        deleted.remove(51);
        assert_eq!(droplet::extract(&deleted).unwrap().trim_start, 86);
        let mut truncated = original.clone();
        truncated.truncate(70);
        assert!(droplet::extract(&truncated).is_err());
        let mut invalid = original.clone();
        invalid[0] = b'N';
        assert!(droplet::extract(&invalid).is_err());
        let mut invalid_umi = original;
        invalid_umi[20] = b'X';
        assert!(droplet::extract(&invalid_umi).is_err());
    }

    #[test]
    fn droplet_correction_rejects_equal_neighbors_but_keeps_exact() {
        let a = "AAAAAAAAAAAAAAAAA".to_string();
        let b = "TGAAAAAAAAAAAAAAA".to_string();
        let index = BarcodeIndex::from_barcodes([a.clone(), b], 17, Modality::Dna);
        assert_eq!(
            index.match_simple(a.as_bytes()),
            BarcodeMatch::Unique { id: &a }
        );
        assert_eq!(
            index.match_simple(b"TAAAAAAAAAAAAAAAA"),
            BarcodeMatch::Ambiguous
        );
        assert_eq!(
            index.match_simple(b"CAAAAAAAAAAAAAAAA"),
            BarcodeMatch::Unique { id: &a }
        );
        assert_eq!(index.match_simple(b"NAAAAAAAAAAAAAAAA"), BarcodeMatch::None);
    }

    #[test]
    fn droplet_fixed_design_keeps_exact_and_rejects_ambiguous_correction() {
        let a = "AAAGAAGAAGAATAGAG".to_string();
        let b = "AAAGAAGAAGAATAGGA".to_string();
        let design = droplet::design_barcodes().unwrap();
        assert_eq!(design.len(), 829_440);
        assert!(design.contains(a.as_bytes()) && design.contains(b.as_bytes()));
        assert!(!design.contains(b"AAAAAAAAAAAAAAAAA".as_slice()));
        let mut index = BarcodeIndex::from_barcodes([a.clone(), b.clone()], 17, Modality::Dna);
        droplet::protect_design(&mut index).unwrap();
        assert_eq!(
            index.match_simple(a.as_bytes()),
            BarcodeMatch::Unique { id: &a }
        );
        assert_eq!(
            index.match_simple(b.as_bytes()),
            BarcodeMatch::Unique { id: &b }
        );
        assert_eq!(
            index.match_simple(b"AAAGAAGAAGAATAGAA"),
            BarcodeMatch::Ambiguous
        );
        assert_eq!(
            index.match_simple(b"CAAGAAGAAGAATAGAG"),
            BarcodeMatch::Unique { id: &a }
        );
        assert_eq!(index.match_simple(b"AAAGAAGAAGAATGAGT"), BarcodeMatch::None);
    }

    use super::*;

    #[test]
    fn input_contract_requires_gzip_fastq_suffix_and_magic_on_both_mates() {
        let root = std::env::temp_dir().join(format!(
            "alopex_gzip_contract.{}.{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir(&root).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut encoder, b"@read\nACGT\n+\nIIII\n").unwrap();
        let gzip = encoder.finish().unwrap();
        let good = root.join("good.fastq.gz");
        let barcode = root.join("Barcode_Map.csv");
        std::fs::write(&good, &gzip).unwrap();
        std::fs::write(&barcode, "DNA_Barcode,PlateID,Cell_Order\nACGTACGT,A1,A1\n").unwrap();
        let make_args = |r1: &Path, r2: &Path| {
            Args::try_parse_from([
                "demux_rs",
                "--r1",
                r1.to_str().unwrap(),
                "--r2",
                r2.to_str().unwrap(),
                "--barcode-map",
                barcode.to_str().unwrap(),
                "--dna-out-dir",
                root.join("dna").to_str().unwrap(),
                "--rna-out-dir",
                root.join("rna").to_str().unwrap(),
                "--json-report",
                root.join("report.json").to_str().unwrap(),
                "--sample-name",
                "AlopexGzip",
                "--mode",
                "dna-only",
            ])
            .unwrap()
        };
        for (name, content, expected) in [
            (
                "plain.fastq",
                b"@read\nACGT\n+\nIIII\n".as_slice(),
                "suffix",
            ),
            ("plain.fasta", b">read\nACGT\n".as_slice(), "suffix"),
            (
                "renamed.fastq.gz",
                b"@read\nACGT\n+\nIIII\n".as_slice(),
                "gzip magic",
            ),
            (
                "renamed_fasta.fq.gz",
                b">read\nACGT\n".as_slice(),
                "gzip magic",
            ),
            ("short.fq.gz", b"\x1f".as_slice(), "gzip header"),
            ("compressed.fasta.gz", gzip.as_slice(), "suffix"),
        ] {
            let bad = root.join(name);
            std::fs::write(&bad, content).unwrap();
            for (r1, r2, label) in [(&bad, &good, "--r1"), (&good, &bad, "--r2")] {
                let error = Config::from_args(make_args(r1, r2))
                    .expect_err("invalid gzip FASTQ")
                    .to_string();
                assert!(
                    error.contains(label) && error.contains(expected),
                    "{name}: {error}"
                );
            }
            assert!(!root.join("report.json").exists());
            assert!(!root.join("dna").exists());
        }
        for name in ["good.fq.gz", "good.FASTQ.GZ"] {
            let path = root.join(name);
            std::fs::write(&path, &gzip).unwrap();
            Config::from_args(make_args(&good, &path)).expect("valid gzip FASTQ");
        }
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn barcode_map_rejects_ambiguous_columns_and_unsafe_plate_ids() {
        let path = std::env::temp_dir().join(format!(
            "dna_pipeline_barcode_audit.{}.csv",
            std::process::id()
        ));
        for content in [
            "DNA_Barcode, dna_barcode ,PlateID,Cell_Order\nACGTACGT,TGCATGCA,A1,1\n",
            "DNA_Barcode,,PlateID,Cell_Order\nACGTACGT,unused,A1,1\n",
            "DNA_Barcode,PlateID,Cell_Order\nACGTACGT,.,1\n",
            "DNA_Barcode,PlateID,Cell_Order\nACGTACGT,..,1\n",
        ] {
            let mut file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
                .expect("isolated barcode fixture");
            std::io::Write::write_all(&mut file, content.as_bytes()).expect("write fixture");
            drop(file);
            let result = load_barcode_map(&path, DemuxMode::DnaOnly);
            std::fs::remove_file(&path).expect("remove fixture");
            assert!(
                result.is_err(),
                "invalid barcode map was accepted: {content}"
            );
        }
    }

    #[test]
    fn sample_name_cannot_escape_the_output_tree() {
        for name in [".", "..", "../outside", "/tmp/outside", "A/B", "A B"] {
            let args = Args::try_parse_from([
                "demux_rs",
                "--r1",
                "r1.fastq.gz",
                "--r2",
                "r2.fastq.gz",
                "--barcode-map",
                "barcodes.csv",
                "--dna-out-dir",
                "dna",
                "--rna-out-dir",
                "rna",
                "--json-report",
                "report.json",
                "--sample-name",
                name,
            ])
            .expect("current CLI arguments");
            let error =
                Config::from_args(args).expect_err("unsafe sample name must fail before I/O");
            assert!(
                error.to_string().contains("--sample-name"),
                "{name}: {error}"
            );
        }
    }

    fn row(dna: &str, rna: &str, order: &str) -> CellRow {
        CellRow {
            dna_barcode: dna.to_string(),
            rna_barcode: rna.to_string(),
            plate_id: format!("P{order}"),
            cell_order: order.to_string(),
        }
    }

    #[test]
    fn paired_read_ids_normalize_mate_suffix_and_description() {
        let normalized = validate_paired_read_ids(b"@read-1/1 lane", b"read-1/2 other")
            .expect("paired IDs should normalize");
        assert_eq!(normalized, b"read-1");
    }

    #[test]
    fn paired_read_id_mismatch_is_rejected() {
        let error = validate_paired_read_ids(b"read-1/1", b"read-2/2")
            .expect_err("different molecule IDs must fail");
        assert!(error.to_string().contains("read-ID mismatch"));
    }

    #[test]
    fn rna_umi_stays_in_the_shared_qname_before_comments() {
        for (left, right) in [
            (&b"read/1"[..], &b"read/2"[..]),
            (&b"read 1:N:0:ATCG"[..], &b"read 2:N:0:ATCG"[..]),
            (&b"read"[..], &b"read"[..]),
        ] {
            for id in [left, right] {
                let mut record = Vec::new();
                write_fastq_record(
                    &mut record,
                    id,
                    b"ACGT",
                    b"IIII",
                    0,
                    Some(b"ACGTACGT"),
                    false,
                )
                .expect("valid RNA FASTQ");
                assert_eq!(
                    record.split(|b| b.is_ascii_whitespace()).next().unwrap(),
                    b"@read:ACGTACGT"
                );
                assert!(record.ends_with(b"\nACGT\n+\nIIII\n"));
            }
        }
    }

    #[test]
    fn protected_barcode_alphabet_and_length_are_enforced() {
        validate_barcode_alphabet("ACGTACGT", "DNA_Barcode", 2).expect("valid protected barcode");
        assert!(validate_barcode_alphabet("ACGTNCGT", "DNA_Barcode", 2).is_err());
        assert!(validate_barcode_alphabet("ACGT", "DNA_Barcode", 2).is_err());
    }

    #[test]
    fn dna_mixed_lengths_are_retained_but_rna_mixed_lengths_fail() {
        let rows = vec![
            row("ACGTACGT", "TGCATGCA", "1"),
            row("ACGTACGTAA", "TGCATGCAAA", "2"),
        ];
        let dna = infer_barcode_lengths(&rows, "DNA", |item| item.dna_barcode.as_str())
            .expect("DNA layout supports registered 8/10bp barcodes");
        assert_eq!(dna, vec![10, 8]);
        assert!(
            infer_uniform_barcode_len(&rows, "RNA", |item| { item.rna_barcode.as_str() }).is_err()
        );
    }

    #[test]
    fn dna_demux_matches_both_registered_8bp_and_10bp_barcodes() {
        init_tables();
        let rows = vec![
            row("ACGTACGT", "TGCATGCA", "1"),
            row("TGCATGCAAA", "ACGTACGTAA", "2"),
        ];
        let lengths = infer_barcode_lengths(&rows, "DNA", |item| item.dna_barcode.as_str())
            .expect("mixed DNA lengths should be supported");
        assert_eq!(lengths, vec![10, 8]);

        let mixed = lengths.len() > 1;
        let indices = lengths
            .iter()
            .map(|&len| {
                BarcodeIndex::from_barcodes(
                    rows.iter()
                        .filter(move |item| !mixed || item.dna_barcode.len() == len)
                        .map(|item| item.dna_barcode.clone()),
                    len,
                    Modality::Dna,
                )
            })
            .collect::<Vec<_>>();
        let rna_index = BarcodeIndex::empty(Modality::Rna);
        let dna_lut = DNA_BISULFITE_TABLE.get().expect("DNA LUT");
        let rna_lut = RNA_BISULFITE_TABLE.get().expect("RNA LUT");

        for barcode in [b"ACGTACGT".as_slice(), b"TGCATGCAAA".as_slice()] {
            let mut read = Vec::new();
            read.extend_from_slice(barcode);
            read.extend_from_slice(ME_SEQ);
            read.extend_from_slice(&vec![b'N'; LEN_GAP_AFTER_ME]);
            read.push(b'G');
            let result = process_read(
                &read,
                &indices,
                &rna_index,
                dna_lut,
                rna_lut,
                0,
                DemuxMode::DnaOnly,
            );
            assert!(result.is_valid, "barcode {:?} should match", barcode);
            assert_eq!(result.fate, TerminalFate::DnaAssigned);
            assert_eq!(result.category, std::str::from_utf8(barcode).unwrap());
        }
    }

    #[test]
    fn taps_me_uses_four_bases_and_keeps_existing_gap_boundary() {
        init_tables();
        let lut = DNA_NATIVE_TABLE.get().unwrap();
        assert_ne!(lut[b'G' as usize], lut[b'A' as usize]);
        let index = BarcodeIndex::from_barcodes(vec!["ACGTACGT".to_string()], 8, Modality::Dna);
        let indices = vec![index];
        let rna = BarcodeIndex::empty(Modality::Rna);
        let mut read = b"ACGTACGT".to_vec();
        read.extend_from_slice(ME_SEQ);
        read.extend_from_slice(b"CCCCCCCCC");
        read.extend_from_slice(b"GACTACGTACGT");
        let result = process_read(
            &read,
            &indices,
            &rna,
            lut,
            RNA_BISULFITE_TABLE.get().unwrap(),
            0,
            DemuxMode::DnaOnlyTaps,
        );
        assert!(result.is_valid);
        assert_eq!(result.fate, TerminalFate::DnaAssigned);
        assert_eq!(result.trim_start, 8 + ME_SEQ.len() + LEN_GAP_AFTER_ME);
        let mut converted = read.clone();
        for base in &mut converted[8..8 + ME_SEQ.len()] {
            if *base == b'G' {
                *base = b'A';
            }
        }
        let result = process_read(
            &converted,
            &indices,
            &rna,
            lut,
            RNA_BISULFITE_TABLE.get().unwrap(),
            0,
            DemuxMode::DnaOnlyTaps,
        );
        assert!(!result.is_valid);
    }

    #[test]
    fn exact_barcode_overrides_ambiguous_one_mismatch_candidate() {
        let index = BarcodeIndex::from_barcodes(
            vec![
                "AAAAAAAA".to_string(),
                "AAAAAAAC".to_string(),
                "AAAAAAAG".to_string(),
            ],
            8,
            Modality::Dna,
        );
        match index.match_simple(b"AAAAAAAG") {
            BarcodeMatch::Unique { id } => {
                assert_eq!(id, "AAAAAAAG");
            }
            other => panic!("expected exact unique match, found {other:?}"),
        }
    }

    #[test]
    fn rna_layout_requires_complete_umi_and_post_umi_spacer() {
        init_tables();
        let index = BarcodeIndex::from_barcodes(vec!["AACGTGAT".to_string()], 8, Modality::Rna);
        let lut = RNA_BISULFITE_TABLE.get().expect("RNA LUT");

        let mut incomplete = Vec::new();
        incomplete.extend_from_slice(TSO_SEQ);
        incomplete.extend_from_slice(b"AACGTGAT");
        incomplete.extend_from_slice(b"ACGTACGT");
        let rejected =
            match_rna_read(&incomplete, &index, lut).expect("layout evidence should exist");
        assert!(!rejected.is_valid);
        assert_eq!(rejected.fate, TerminalFate::Unmatched);

        let mut complete = incomplete.clone();
        complete.extend_from_slice(b"NNNNN");
        complete.extend_from_slice(b"GATTACA");
        let assigned = match_rna_read(&complete, &index, lut).expect("complete RNA layout");
        assert!(assigned.is_valid);
        assert_eq!(assigned.fate, TerminalFate::RnaAssigned);
        assert_eq!(
            assigned.trim_start,
            TSO_SEQ.len() + 8 + LEN_UMI + LEN_SPACER
        );
    }

    #[test]
    fn observed_n_is_not_imputed_by_one_substitution_matching() {
        let index = BarcodeIndex::from_barcodes(vec!["ACGTACGT".to_string()], 8, Modality::Dna);
        assert_eq!(index.match_simple(b"ACGTACGN"), BarcodeMatch::None);
    }

    #[test]
    fn rna_layout_requires_at_least_one_payload_base() {
        init_tables();
        let index = BarcodeIndex::from_barcodes(vec!["AACGTGAT".to_string()], 8, Modality::Rna);
        let lut = RNA_BISULFITE_TABLE.get().expect("RNA LUT");

        let mut no_payload = Vec::new();
        no_payload.extend_from_slice(TSO_SEQ);
        no_payload.extend_from_slice(b"AACGTGAT");
        no_payload.extend_from_slice(b"ACGTACGT");
        no_payload.extend_from_slice(b"NNNNN");
        let rejected =
            match_rna_read(&no_payload, &index, lut).expect("layout evidence should exist");
        assert!(!rejected.is_valid);
        assert_eq!(rejected.fate, TerminalFate::Unmatched);

        no_payload.push(b'G');
        let assigned =
            match_rna_read(&no_payload, &index, lut).expect("payload should complete layout");
        assert!(assigned.is_valid);
        assert_eq!(assigned.fate, TerminalFate::RnaAssigned);
    }

    #[test]
    fn dna_layout_requires_at_least_one_insert_base_after_gap() {
        init_tables();
        let index = BarcodeIndex::from_barcodes(vec!["ACGTACGT".to_string()], 8, Modality::Dna);
        let lut = DNA_BISULFITE_TABLE.get().expect("DNA LUT");

        let mut no_insert = Vec::new();
        no_insert.extend_from_slice(b"ACGTACGT");
        no_insert.extend_from_slice(ME_SEQ);
        no_insert.extend_from_slice(&vec![b'N'; LEN_GAP_AFTER_ME]);
        let rejected =
            match_dna_read(&no_insert, &index, lut, 0).expect("layout evidence should exist");
        assert!(!rejected.is_valid);
        assert_eq!(rejected.fate, TerminalFate::Unmatched);

        no_insert.push(b'G');
        let assigned =
            match_dna_read(&no_insert, &index, lut, 0).expect("insert should complete layout");
        assert!(assigned.is_valid);
        assert_eq!(assigned.fate, TerminalFate::DnaAssigned);
    }
}

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
        if let Some(existing) = order_map.insert(row.cell_order.clone(), row.dna_barcode.clone()) {
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

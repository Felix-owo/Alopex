use super::*;

pub(super) fn process_read<'a>(
    r1: &'a [u8],
    dna_indices: &'a [BarcodeIndex],
    dna_lut: &'static [u8; 256],
    dna_w_spacer_len: usize,
) -> ProcResult<'a> {
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
                return Some(ProcResult::assigned(id, payload_start));
            }
        }
        retain_best_rejection(
            &mut rejection,
            rejected_barcode_attempt(barcode_match),
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
                let payload_start = alignment.end + LEN_GAP_AFTER_ME;
                if payload_start < r1.len() {
                    if let BarcodeMatch::Unique { id } = barcode_match {
                        return Some(ProcResult::assigned(id, payload_start));
                    }
                }
                retain_best_rejection(
                    &mut rejection,
                    rejected_barcode_attempt(barcode_match),
                );
            }
        }
    }

    rejection
}

pub(super) fn rejected_barcode_attempt<'a>(
    barcode_match: BarcodeMatch<'_>,
) -> ProcResult<'a> {
    let fate = match barcode_match {
        BarcodeMatch::Ambiguous => TerminalFate::Ambiguous,
        BarcodeMatch::Unique { .. } | BarcodeMatch::None => TerminalFate::Unmatched,
    };
    ProcResult::rejected(fate)
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

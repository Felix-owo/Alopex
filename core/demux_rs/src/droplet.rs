pub(crate) mod evidence;
pub(crate) mod replay;
mod spool;
pub(crate) mod stage;

use super::*;
use std::io::{BufRead, BufReader};

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
    if !barcode.iter().all(|b| b"ACGT".contains(b)) || !umi.iter().all(|b| b"ACGTN".contains(b)) {
        return Err(TerminalFate::Unmatched);
    }
    Ok(Layout {
        barcode,
        umi,
        barcode_start: start - 29,
        trim_start: end + 9,
    })
}

#[derive(Default, Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct CountMetrics {
    pub(super) input_reads: u64,
    pub(super) structured_reads: u64,
    pub(super) ambiguous: u64,
    pub(super) unmatched: u64,
    pub(super) short: u64,
    pub(super) barcode_with_c: u64,
    pub(super) umi_with_n: u64,
    pub(super) unique_barcodes: usize,
}

#[derive(Default)]
pub(super) struct BarcodeCount {
    pub(super) count: u64,
    pub(super) below_q20: [u64; 17],
    pub(super) q30: [u64; 17],
}

pub(super) fn write_count_row<W: std::io::Write>(
    out: &mut W,
    barcode: &[u8],
    count: &BarcodeCount,
) -> Result<()> {
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
        std::str::from_utf8(barcode)?,
        count.count,
        list(&count.below_q20),
        list(&count.q30)
    )?;
    Ok(())
}

const MOLECULE_CAP: usize = 4096;

pub(super) fn molecule_key(layout: &Layout<'_>, seq: &[u8]) -> Option<u128> {
    let insert = seq.get(layout.trim_start..layout.trim_start + 32)?;
    encode_molecule(layout.umi, insert)
}

pub(super) fn encode_molecule(umi: &[u8], insert: &[u8]) -> Option<u128> {
    umi.iter().chain(insert).try_fold(0u128, |key, base| {
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

pub(super) fn molecule_rank(key: u128) -> (u64, u128) {
    let mix = |mut x: u64| {
        x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
        x ^ (x >> 31)
    };
    (mix((key as u64) ^ mix((key >> 64) as u64)), key)
}

#[derive(Default)]
pub(super) struct Molecules {
    pub(super) reads: u64,
    pub(super) saturated: bool,
    pub(super) keys: std::collections::BTreeSet<(u64, u128)>,
}

impl Molecules {
    pub(super) fn add(&mut self, key: u128) {
        self.reads += 1;
        let ranked = molecule_rank(key);
        if self.keys.len() < MOLECULE_CAP || self.keys.last().is_some_and(|last| ranked <= *last) {
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

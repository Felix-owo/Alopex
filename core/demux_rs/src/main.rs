//! 解复用入口、共享记录类型与孔板处理；Droplet阶段实现由职责模块承载。
mod barcode_map;
mod cli;
mod compression;
mod droplet;
mod fastq_io;
mod matcher;
mod output;
mod writeback;

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
            ensure_fastq_input(path, "--r1")?;
        }
        for path in &args.r2 {
            ensure_fastq_input(path, "--r2")?;
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

    if std::env::args().nth(1).as_deref() == Some("droplet-stage") {
        return droplet::stage::run();
    }
    let args = Args::parse();
    let droplet_spool = args.droplet_spool.clone();
    anyhow::ensure!(
        droplet_spool.is_none() || args.mode == DemuxMode::DnaOnlyDroplet,
        "staged reads require dna-only-droplet"
    );
    anyhow::ensure!(
        args.mode != DemuxMode::DnaOnlyDroplet || droplet_spool.is_some(),
        "dna-only-droplet requires --droplet-spool from droplet-stage"
    );
    if args.threads == 0 {
        anyhow::bail!("--threads must be greater than 0");
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(
            args.threads
                .saturating_sub(if droplet_spool.is_some() { 4 } else { 0 })
                .max(1),
        )
        .build_global()
        .context("Cannot initialize Rayon global thread pool")?;
    if args.evidence_pairs.is_some() {
        return droplet::evidence::run(&args, droplet_spool.as_deref().unwrap());
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

    if let Some(spool) = droplet_spool {
        return droplet::replay::run(&config, cell_rows, dna_indices, paths, &spool);
    }
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

    let (tx_completion, rx_completion) = bounded(1);
    let writer_handle = spawn_writer(rx_write, tx_acc_recycle, paths, cell_rows, rx_completion);

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
    let _ = tx_completion.send(read_result.is_ok() && worker_result.is_ok());
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
        let result = process_read(
            r1_seq,
            dna_indices,
            rna_idx,
            dna_lut,
            rna_lut,
            dna_w_spacer_len,
            demux_mode,
        );
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
            &input.r1_ids[r1_id_start..r1_id_end],
            r1_seq,
            &input.r1_quals[r1_start..r1_end],
            result.trim_start,
            result.output_umi(),
            result.is_valid,
            true,
        )?;
        write_fastq_record(
            &mut bucket.r2_raw,
            &input.r2_ids[r2_id_start..r2_id_end],
            &input.r2_seqs[r2_start..r2_end],
            &input.r2_quals[r2_start..r2_end],
            0,
            result.output_umi(),
            false,
            true,
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

fn ensure_fastq_input(path: &Path, label: &str) -> Result<()> {
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
        let normalized = validate_paired_read_ids(b"read-1/1 lane", b"read-1/2 other")
            .expect("paired IDs should normalize");
        assert_eq!(normalized, b"read-1");
    }

    #[test]
    fn parsed_fastq_ids_keep_literal_at_and_reject_blank_identifiers() {
        let input = b"@read/1 lane\nACGT\n+\nIIII\n@@read/1 lane\nACGT\n+\nIIII\n@read/1/1 lane\nACGT\n+\nIIII\n";
        let mut reader = needletail::parse_fastx_reader(&input[..]).unwrap();
        let mut names = Vec::new();
        while let Some(record) = reader.next() {
            names.push(normalize_read_id(record.unwrap().id()).unwrap().to_vec());
        }
        assert_eq!(names, [b"read".to_vec(), b"@read".to_vec(), b"read/1".to_vec()]);
        assert!(validate_paired_read_ids(b"@read/1", b"read/2").is_err());
        assert!(validate_paired_read_ids(b"read/1", b"@read/2").is_err());
        for id in [b"".as_slice(), b" read", b"\tread", b"/1"] {
            assert!(normalize_read_id(id).is_err());
        }
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
                    true,
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
    fn dna_indel_fallback_trims_at_actual_anchor_end() {
        init_tables();
        let index = BarcodeIndex::from_barcodes(["ACGTACGT".to_string()], 8, Modality::Dna);
        for lut in [
            DNA_NATIVE_TABLE.get().unwrap(),
            DNA_BISULFITE_TABLE.get().unwrap(),
        ] {
            for insertion in [true, false] {
                let mut anchor = ME_SEQ.to_vec();
                if insertion {
                    anchor.insert(8, b'C');
                } else {
                    anchor.remove(8);
                }
                let read = [b"ACGTACGT".as_slice(), &anchor, b"NNNNNNNNN", b"GACTACGT"].concat();
                let assigned = match_dna_read(&read, &index, lut, 0).unwrap();
                assert!(assigned.is_valid);
                assert_eq!(
                    &read[assigned.trim_start..],
                    b"GACTACGT",
                    "insertion={insertion}"
                );
            }
        }
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

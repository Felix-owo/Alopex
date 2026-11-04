use clap::{Parser, ValueEnum};
use std::path::PathBuf;

#[derive(Parser, Debug, Clone)]
#[command(author, version, about = "DNA/RNA FASTQ 高通量解复用器")]
pub(crate) struct Args {
    #[arg(
        long,
        required_unless_present = "evidence_pairs",
        help = "R1 FASTQ；可按 lane/chunk 顺序重复传入。"
    )]
    pub(crate) r1: Vec<PathBuf>,

    #[arg(
        long,
        required_unless_present = "evidence_pairs",
        help = "R2 FASTQ；数量和顺序必须与 --r1 一一对应。"
    )]
    pub(crate) r2: Vec<PathBuf>,

    #[arg(long, help = "Droplet 配对暂存目录；该模式的重放与分子证据必需。")]
    pub(crate) droplet_spool: Option<PathBuf>,

    #[arg(long, requires_all = ["evidence_output", "droplet_spool"], conflicts_with_all = ["r1", "r2"])]
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

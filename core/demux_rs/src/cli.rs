use clap::{Parser, ValueEnum};
use std::path::{Component, Path, PathBuf};
use std::os::unix::fs::MetadataExt;
use crate::{Context, Result};

fn resolve_output_path(path: &Path) -> Result<PathBuf> {
    let mut resolved = std::env::current_dir()?;
    let mut components = path.components().peekable();
    while let Some(component) = components.next() {
        match component {
            Component::RootDir => resolved = PathBuf::from("/"),
            Component::CurDir => {},
            Component::ParentDir => { resolved.pop(); },
            Component::Normal(name) => {
                resolved.push(name);
                match std::fs::symlink_metadata(&resolved) {
                    Ok(metadata) if metadata.file_type().is_symlink() => {
                        anyhow::ensure!(components.peek().is_some(),
                            "Output must not be a symlink: {}", path.display());
                        resolved = resolved.canonicalize().with_context(||
                            format!("Cannot resolve output symlink: {}", resolved.display()))?;
                        anyhow::ensure!(std::fs::metadata(&resolved)?.is_dir(),
                            "Output ancestor must be a directory: {}", resolved.display());
                    },
                    Ok(metadata) => {
                        anyhow::ensure!(components.peek().is_none() || metadata.is_dir(),
                            "Output ancestor must be a directory: {}", resolved.display());
                    },
                    Err(error) if error.kind() == std::io::ErrorKind::NotFound => {},
                    Err(error) => return Err(error).with_context(||
                        format!("Cannot inspect output path: {}", resolved.display())),
                }
            },
            Component::Prefix(_) => anyhow::bail!("Unsupported output path: {}", path.display()),
        }
    }
    Ok(resolved)
}

/// 在任何输出创建、清理或截断之前验证整组路径；主流程、stage和evidence共用。
pub(crate) fn validate_io_paths(
    inputs: &[(&str, &Path)], directories: &[(&str, &Path)], files: &[(&str, &Path)],
) -> Result<()> {
    let protected = inputs.iter().map(|(label, path)| {
        let resolved = path.canonicalize().with_context(||
            format!("Cannot resolve {label}: {}", path.display()))?;
        let directory = std::fs::metadata(&resolved)?.is_dir();
        Ok((*label, resolved, directory))
    }).collect::<Result<Vec<_>>>()?;
    let mut outputs: Vec<(&str, PathBuf, bool)> = Vec::new();
    for (directory, paths) in [(true, directories), (false, files)] {
        for (label, path) in paths {
            let resolved = resolve_output_path(path)?;
            match std::fs::symlink_metadata(&resolved) {
                Ok(metadata) => {
                    anyhow::ensure!(!metadata.file_type().is_symlink(),
                        "{label} must not be a symlink: {}", path.display());
                    anyhow::ensure!(if directory { metadata.is_dir() } else { metadata.is_file() },
                        "{label} has an invalid output type: {}", path.display());
                    // Pipeline的report可预先由mktemp创建；不得沿共享inode截断输入或备份。
                    anyhow::ensure!(directory || metadata.nlink() == 1,
                        "{label} must not be hard-linked: {}", path.display());
                },
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {},
                Err(error) => return Err(error).with_context(||
                    format!("Cannot inspect {label}: {}", path.display())),
            }
            for (input_label, input, input_directory) in &protected {
                anyhow::ensure!(resolved != *input
                    && !(directory && input.starts_with(&resolved))
                    && !(*input_directory && resolved.starts_with(input)),
                    "{label} overlaps protected {input_label}: {}", path.display());
            }
            for (other_label, other, other_directory) in &outputs {
                anyhow::ensure!(resolved != *other
                    && !(directory && other.starts_with(&resolved))
                    && !(*other_directory && resolved.starts_with(other)),
                    "{label} overlaps {other_label}: {}", path.display());
            }
            outputs.push((*label, resolved, directory));
        }
    }
    Ok(())
}

#[test]
fn io_plan_rejects_aliases_without_creating_or_changing_files() {
    use std::os::unix::fs::symlink;
    let root = std::env::temp_dir().join(format!("alopex_io_plan.{}.{}", std::process::id(),
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
    let raw = root.join("raw");
    std::fs::create_dir_all(&raw).unwrap();
    let input = raw.join("R1.fastq.gz");
    std::fs::write(&input, b"original molecule records").unwrap();
    let redirect = root.join("redirect");
    symlink(&raw, &redirect).unwrap();
    let missing = root.join("not-created");
    let directory = missing.join("../redirect");
    let report = root.join("report.json");
    let inputs = [("input", input.as_path())];
    assert!(validate_io_paths(&inputs, &[("DNA", &directory)], &[]).is_err());
    assert!(validate_io_paths(&inputs, &[], &[("report", &input)]).is_err());
    let linked = root.join("linked-report.json");
    symlink(&input, &linked).unwrap();
    assert!(validate_io_paths(&inputs, &[], &[("report", &linked)]).is_err());
    let hard = root.join("shared-report.json");
    std::fs::hard_link(&input, &hard).unwrap();
    assert!(validate_io_paths(&inputs, &[], &[("report", &hard)]).is_err());
    assert!(validate_io_paths(&inputs, &[], &[("report", &missing.join("../shared-report.json"))]).is_err());
    let not_directory = root.join("file");
    std::fs::write(&not_directory, b"not a directory").unwrap();
    assert!(validate_io_paths(&inputs, &[("RNA", &missing.join("../file/../RNA"))], &[]).is_err());
    let output = missing.join("DNA");
    assert!(validate_io_paths(&inputs, &[("DNA", &output)],
        &[("report", &output.join("report.json"))]).is_err());
    assert!(validate_io_paths(&inputs, &[], &[("counts", &report), ("metrics", &report)]).is_err());
    assert!(validate_io_paths(&[("spool", raw.as_path())], &[], &[("report", &raw.join("manifest.json"))]).is_err());
    // Pipeline独占mktemp report和正常独立目录仍合法；校验不提前创建目录或截断文件。
    std::fs::write(&report, b"").unwrap();
    validate_io_paths(&inputs, &[("DNA", &output)], &[("report", &report)]).unwrap();
    assert!(!missing.exists());
    assert_eq!(std::fs::read(&input).unwrap(), b"original molecule records");
    assert!(std::fs::read(&report).unwrap().is_empty());
    std::fs::remove_dir_all(root).unwrap();
}

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

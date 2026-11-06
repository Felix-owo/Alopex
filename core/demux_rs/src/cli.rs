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

struct ProtectedInput<'a> {
    label: &'a str,
    path: PathBuf,
    subtree: bool,
    identity: (u64, u64),
}

struct PlannedOutput<'a> {
    label: &'a str,
    path: PathBuf,
    directory: bool,
    identity: Option<(u64, u64)>,
    parents: Vec<(u64, u64)>,
}

/// 在目标文件系统独占临时目录中模拟尚未存在的命名空间；不创建最终输出。
pub(crate) fn validate_prospective_paths(
    directories: &[&Path], files: &[&Path], replaced_root: Option<&Path>,
) -> Result<()> {
    use std::collections::BTreeMap;
    use std::fs::{self, OpenOptions};
    use std::os::unix::fs::DirBuilderExt;
    let replacement = replaced_root.map(|root| resolve_output_path(root).map(|path| (root, path))).transpose()?;
    let mut groups: BTreeMap<(u64, u64), (PathBuf, Vec<(PathBuf, bool)>)> = BTreeMap::new();
    for (directory, paths) in [(true, directories), (false, files)] {
        for path in paths {
            // 独占根将被删除重建；旧根权限/目录标志/子链接不定义新根的名称语义。
            let path = if let Some((root, resolved)) = &replacement {
                resolved.join(path.strip_prefix(root)?)
            } else { resolve_output_path(path)? };
            let namespace = replacement.as_ref().map(|(_, root)| root.as_path()).unwrap_or(&path);
            let mut parent = namespace.parent().context("Output has no parent")?;
            while !parent.exists() {
                parent = parent.parent().context("Output has no existing ancestor")?;
            }
            let metadata = fs::metadata(parent)?;
            anyhow::ensure!(metadata.is_dir(), "Output ancestor must be a directory: {}", parent.display());
            let suffix = path.strip_prefix(parent)?.to_path_buf();
            groups.entry((metadata.dev(), metadata.ino()))
                .or_insert_with(|| (parent.to_path_buf(), Vec::new())).1.push((suffix, directory));
        }
    }
    for (parent_identity, (parent, members)) in groups {
        if members.len() < 2 { continue; }
        let parent = parent.canonicalize()?;
        let nonce = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)?.as_nanos();
        let probe = (0..128).find_map(|attempt| {
            let path = parent.join(format!(".alopex_output_names.{}.{nonce}.{attempt}", std::process::id()));
            match fs::DirBuilder::new().mode(0o700).create(&path) {
                Ok(()) => Some(Ok(path)),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => None,
                Err(error) => Some(Err(error)),
            }
        }).context("Cannot allocate exclusive output namespace probe")??;
        let metadata = fs::symlink_metadata(&probe)?;
        let probe_identity = (metadata.dev(), metadata.ino());
        let result = (|| -> Result<()> {
            let mut outputs = crate::AHashSet::new();
            let mut output_directories = crate::AHashSet::new();
            let mut output_ancestors = crate::AHashSet::new();
            for (suffix, directory) in &members {
                let path = probe.join(suffix);
                fs::create_dir_all(path.parent().context("Probe output has no parent")?)?;
                if *directory { fs::create_dir_all(&path)?; }
                else { OpenOptions::new().write(true).create_new(true).open(&path)
                    .context("Output names overlap on the target filesystem")?; }
                let metadata = fs::metadata(&path)?;
                let identity = (metadata.dev(), metadata.ino());
                let parents: Vec<_> = path.ancestors().skip(1).map(|p| {
                    fs::metadata(p).map(|m| (m.dev(), m.ino()))
                }).collect::<std::io::Result<_>>()?;
                anyhow::ensure!(!outputs.contains(&identity)
                    && !(*directory && output_ancestors.contains(&identity))
                    && !parents.iter().any(|id| output_directories.contains(id)),
                    "Output names overlap on the target filesystem: {}", suffix.display());
                outputs.insert(identity);
                if *directory { output_directories.insert(identity); }
                output_ancestors.extend(parents);
            }
            Ok(())
        })();
        let cleanup = (|| -> Result<()> {
            // 清理只认本次父目录和临时根身份，不按前缀搜索，也不跟随替换后的根。
            let parent_metadata = fs::symlink_metadata(&parent)?;
            let metadata = fs::symlink_metadata(&probe)?;
            anyhow::ensure!(!parent_metadata.file_type().is_symlink() && !metadata.file_type().is_symlink()
                && (parent_metadata.dev(), parent_metadata.ino()) == parent_identity
                && (metadata.dev(), metadata.ino()) == probe_identity,
                "Output namespace probe changed before cleanup: {}", probe.display());
            fs::remove_dir_all(&probe)?;
            Ok(())
        })();
        if let Err(error) = cleanup {
            if result.is_ok() { return Err(error); }
            log::warn!("Output namespace probe cleanup failed; preserved {}: {error:#}", probe.display());
        }
        result?;
    }
    Ok(())
}

fn protect_input_entries<'a>(
    path: &Path, base: &Path, label: &'a str, protected: &mut Vec<ProtectedInput<'a>>,
) -> Result<PathBuf> {
    let mut resolved = base.to_path_buf();
    let mut components = path.components().peekable();
    while let Some(component) = components.next() {
        match component {
            Component::RootDir => resolved = PathBuf::from("/"),
            Component::CurDir => {},
            Component::ParentDir => { resolved.pop(); },
            Component::Normal(name) => {
                resolved.push(name);
                // 目录重置不能删掉输入所经的链接；只保canonical目标会漏掉入口及中间跳。
                let mut metadata = std::fs::symlink_metadata(&resolved)?;
                protected.push(ProtectedInput { label, path: resolved.clone(), subtree: false,
                    identity: (metadata.dev(), metadata.ino()) });
                if metadata.file_type().is_symlink() {
                    // 沿用文件系统的循环/无效链接拒绝，再逐跳记录实际读取路径。
                    resolved.canonicalize()?;
                    let target = std::fs::read_link(&resolved)?;
                    let parent = resolved.parent().context("Input symlink has no parent")?;
                    resolved = protect_input_entries(&target, parent, label, protected)?;
                    metadata = std::fs::metadata(&resolved)?;
                }
                anyhow::ensure!(components.peek().is_none() || metadata.is_dir(),
                    "{label} input ancestor must be a directory: {}", path.display());
            },
            Component::Prefix(_) => anyhow::bail!("Unsupported input path: {}", path.display()),
        }
    }
    Ok(resolved)
}

/// 在任何输出创建、清理或截断之前验证整组路径；主流程、stage和evidence共用。
pub(crate) fn validate_io_paths(
    inputs: &[(&str, &Path)], directories: &[(&str, &Path)], files: &[(&str, &Path)],
) -> Result<()> {
    let mut protected = Vec::new();
    for (label, path) in inputs {
        let resolved = path.canonicalize().with_context(||
            format!("Cannot resolve {label}: {}", path.display()))?;
        let metadata = std::fs::metadata(&resolved)?;
        let traced = protect_input_entries(path, &std::env::current_dir()?, label, &mut protected)?;
        let traced_metadata = std::fs::metadata(&traced)?;
        anyhow::ensure!((metadata.dev(), metadata.ino()) == (traced_metadata.dev(), traced_metadata.ino()),
            "{label} changed while resolving input: {}", path.display());
        protected.push(ProtectedInput { label, path: resolved, subtree: metadata.is_dir(),
            identity: (metadata.dev(), metadata.ino()) });
    }
    let mut outputs: Vec<PlannedOutput<'_>> = Vec::new();
    for (directory, paths) in [(true, directories), (false, files)] {
        for (label, path) in paths {
            let resolved = resolve_output_path(path)?;
            let mut identity = None;
            match std::fs::symlink_metadata(&resolved) {
                Ok(metadata) => {
                    identity = Some((metadata.dev(), metadata.ino()));
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
            let mut parents = Vec::new();
            for parent in resolved.ancestors().skip(1) {
                match std::fs::metadata(parent) {
                    Ok(metadata) => {
                        anyhow::ensure!(metadata.is_dir(), "Output ancestor must be a directory: {}", parent.display());
                        parents.push((metadata.dev(), metadata.ino()));
                    },
                    Err(error) if error.kind() == std::io::ErrorKind::NotFound => {},
                    Err(error) => return Err(error).with_context(|| format!("Cannot inspect output ancestor: {}", parent.display())),
                }
            }
            for input in &protected {
                anyhow::ensure!(resolved != input.path
                    && !(directory && input.path.starts_with(&resolved))
                    && !(input.subtree && resolved.starts_with(&input.path))
                    && identity != Some(input.identity)
                    && !(input.subtree && parents.contains(&input.identity)),
                    "{label} overlaps protected {}: {}", input.label, path.display());
            }
            for other in &outputs {
                anyhow::ensure!(resolved != other.path
                    && !(directory && other.path.starts_with(&resolved))
                    && !(other.directory && resolved.starts_with(&other.path))
                    && !(identity.is_some() && identity == other.identity)
                    && !(directory && identity.is_some_and(|id| other.parents.contains(&id)))
                    && !(other.directory && other.identity.is_some_and(|id| parents.contains(&id))),
                    "{label} overlaps {}: {}", other.label, path.display());
            }
            outputs.push(PlannedOutput { label, path: resolved, directory, identity, parents });
        }
    }
    validate_prospective_paths(
        &outputs.iter().filter(|p| p.directory).map(|p| p.path.as_path()).collect::<Vec<_>>(),
        &outputs.iter().filter(|p| !p.directory).map(|p| p.path.as_path()).collect::<Vec<_>>(),
        None,
    )
}

#[test]
fn input_entries_and_existing_aliases_are_protected_without_mutation() {
    use std::os::unix::fs::symlink;
    let root = std::env::temp_dir().join(format!("alopex_input_entries.{}.{}", std::process::id(),
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
    let raw = root.join("Raw");
    let output = root.join("Output");
    std::fs::create_dir_all(&raw).unwrap();
    std::fs::create_dir(&output).unwrap();
    let input = raw.join("R1.fastq.gz");
    std::fs::write(&input, b"original input").unwrap();
    let marker = output.join("old");
    std::fs::write(&marker, b"old output").unwrap();
    let middle = output.join("middle");
    symlink(&input, &middle).unwrap();
    let entry = root.join("entry");
    symlink("Output/middle", &entry).unwrap();
    assert!(validate_io_paths(&[("raw", &entry)], &[("DNA", &output)], &[]).is_err());
    let middle_directory = output.join("reads");
    symlink(&raw, &middle_directory).unwrap();
    let directory_entry = root.join("reads");
    symlink("Output/reads", &directory_entry).unwrap();
    assert!(validate_io_paths(&[("raw", &directory_entry.join("R1.fastq.gz"))], &[("DNA", &output)], &[]).is_err());
    // 合法相对链接链及重复经过同一目录链接仍可读，不把再次出现的链接误判为环。
    symlink("Raw", root.join("legal")).unwrap();
    let legal = root.join("legal/../legal/R1.fastq.gz");
    let report = root.join("report.json");
    validate_io_paths(&[("raw", &legal)], &[("DNA", &root.join("new"))], &[("report", &report)]).unwrap();
    symlink("cycle", root.join("cycle")).unwrap();
    assert!(validate_io_paths(&[("raw", &root.join("cycle"))], &[], &[("report", &report)]).is_err());
    // 大小写别名仅在当前文件系统确实支持时检查；其他平台仍执行上面的完整路径用例。
    let raw_alias = root.join("raw");
    if raw_alias.exists() {
        assert!(validate_io_paths(&[("raw", &input)], &[("DNA", &raw_alias)], &[]).is_err());
        assert!(validate_io_paths(&[("spool", &raw)], &[], &[("evidence", &raw_alias.join("new.json"))]).is_err());
        assert!(validate_io_paths(&[("raw", &input)], &[("DNA", &output), ("RNA", &root.join("output"))], &[]).is_err());
        assert!(validate_io_paths(&[("raw", &input)], &[("DNA", &output)], &[("report", &root.join("output/report.json"))]).is_err());
    }
    assert_eq!(std::fs::read(&input).unwrap(), b"original input");
    assert_eq!(std::fs::read(&marker).unwrap(), b"old output");
    assert_eq!(std::fs::read_link(&entry).unwrap(), PathBuf::from("Output/middle"));
    assert!(!root.join("new").exists());
    assert!(!report.exists());
    std::fs::remove_dir_all(root).unwrap();
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

//! 完整双端暂存、结构计数和manifest最后提交；不读取最终called名单。
use super::spool::*;
use crate::compression::compress_stage_shards;
use crate::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

#[derive(Parser, Debug)]
struct StageArgs {
    #[arg(long, required = true)]
    r1: Vec<PathBuf>,
    #[arg(long, required = true)]
    r2: Vec<PathBuf>,
    #[arg(long)]
    spool: PathBuf,
    #[arg(long)]
    counts_output: PathBuf,
    #[arg(long)]
    count_metrics: PathBuf,
    #[arg(long, default_value_t = 4)]
    threads: usize,
    #[arg(long, default_value_t = 64)]
    shards: usize,
    #[arg(long, default_value_t = 100_000)]
    batch_pairs: usize,
    #[arg(long)]
    decoder_python: PathBuf,
}
enum QualityCounts {
    Small([u8; 34]),
    Wide(Box<[u64; 34]>),
}

struct CompactCount {
    count: u64,
    quality: QualityCounts,
}

impl Default for CompactCount {
    fn default() -> Self {
        Self {
            count: 0,
            quality: QualityCounts::Small([0; 34]),
        }
    }
}

impl CompactCount {
    fn merge(&mut self, added: &droplet::BarcodeCount) -> Result<()> {
        let total = self
            .count
            .checked_add(added.count)
            .context("barcode count overflow")?;
        if let QualityCounts::Small(values) = &self.quality {
            let increments = added.below_q20.iter().chain(added.q30.iter());
            if values
                .iter()
                .zip(increments)
                .any(|(old, add)| *add > u64::from(u8::MAX - *old))
            {
                self.quality =
                    QualityCounts::Wide(Box::new(std::array::from_fn(|i| u64::from(values[i]))));
            }
        }
        match &mut self.quality {
            QualityCounts::Small(values) => {
                for (old, add) in values
                    .iter_mut()
                    .zip(added.below_q20.iter().chain(added.q30.iter()))
                {
                    *old += u8::try_from(*add)?;
                }
            }
            QualityCounts::Wide(values) => {
                for (old, add) in values
                    .iter_mut()
                    .zip(added.below_q20.iter().chain(added.q30.iter()))
                {
                    *old = old
                        .checked_add(*add)
                        .context("barcode quality count overflow")?;
                }
            }
        }
        self.count = total;
        Ok(())
    }

    fn expanded(&self) -> droplet::BarcodeCount {
        let get = |i| match &self.quality {
            QualityCounts::Small(values) => u64::from(values[i]),
            QualityCounts::Wide(values) => values[i],
        };
        droplet::BarcodeCount {
            count: self.count,
            below_q20: std::array::from_fn(|i| get(i)),
            q30: std::array::from_fn(|i| get(i + 17)),
        }
    }
}

struct CountTable {
    parts: [AHashMap<[u8; 17], CompactCount>; 64],
}

impl Default for CountTable {
    fn default() -> Self {
        Self {
            parts: std::array::from_fn(|_| AHashMap::new()),
        }
    }
}

impl CountTable {
    fn merge(&mut self, barcode: [u8; 17], value: &droplet::BarcodeCount) -> Result<()> {
        self.parts[shard(&barcode, self.parts.len())]
            .entry(barcode)
            .or_default()
            .merge(value)
    }
    fn get(&self, barcode: &[u8; 17]) -> &CompactCount {
        &self.parts[shard(barcode, self.parts.len())][barcode]
    }
    fn len(&self) -> usize {
        self.parts.iter().map(|part| part.len()).sum()
    }
    fn iter(&self) -> impl Iterator<Item = (&[u8; 17], &CompactCount)> {
        self.parts.iter().flat_map(|part| part.iter())
    }
}

fn write_compact_counts(
    counts: CountTable,
    mut metrics: droplet::CountMetrics,
    counts_path: &Path,
    metrics_path: &Path,
) -> Result<()> {
    metrics.unique_barcodes = counts.len();
    let mut rows = counts
        .iter()
        .map(|(key, count)| (*key, count.count))
        .collect::<Vec<_>>();
    rows.par_sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let mut out = BufWriter::with_capacity(1024 * 1024, File::create(counts_path)?);
    writeln!(out, "barcode\tcount\tbelow_q20\tq30")?;
    for (barcode, _) in rows {
        droplet::write_count_row(&mut out, &barcode, &counts.get(&barcode).expanded())?;
    }
    out.flush()?;
    std::fs::write(metrics_path, serde_json::to_vec_pretty(&metrics)?)?;
    Ok(())
}

#[test]
fn compact_count_promotes_losslessly_and_rejects_overflow() {
    assert!(std::mem::size_of::<CompactCount>() <= 48);
    let mut value = CompactCount::default();
    value
        .merge(&droplet::BarcodeCount {
            count: 255,
            below_q20: [254; 17],
            q30: [1; 17],
        })
        .unwrap();
    assert!(matches!(value.quality, QualityCounts::Small(_)));
    value
        .merge(&droplet::BarcodeCount {
            count: 5,
            below_q20: [2; 17],
            q30: [3; 17],
        })
        .unwrap();
    assert!(matches!(value.quality, QualityCounts::Wide(_)));
    let expanded = value.expanded();
    assert_eq!(expanded.count, 260);
    assert_eq!(expanded.below_q20, [256; 17]);
    assert_eq!(expanded.q30, [4; 17]);
    value.count = u64::MAX;
    assert!(value
        .merge(&droplet::BarcodeCount {
            count: 1,
            ..Default::default()
        })
        .is_err());
}

#[test]
fn compact_count_output_preserves_u64_quality_counts_and_tie_order() {
    let root = std::env::temp_dir().join(format!("demux-compact-counts-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let mut counts = CountTable::default();
    let values = [1u64, 254, 255, 256, 999, 999, u32::MAX as u64 + 10];
    for (index, count) in values.into_iter().enumerate() {
        let mut key = [b'A'; 17];
        key[16] = b'A' + index as u8;
        counts
            .merge(
                key,
                &droplet::BarcodeCount {
                    count,
                    below_q20: [count / 3; 17],
                    q30: [count - count / 3; 17],
                },
            )
            .unwrap();
    }
    write_compact_counts(
        counts,
        Default::default(),
        &root.join("counts.tsv"),
        &root.join("metrics.json"),
    )
    .unwrap();
    let output = std::fs::read_to_string(root.join("counts.tsv")).unwrap();
    let mut lines = output.lines();
    assert_eq!(lines.next(), Some("barcode\tcount\tbelow_q20\tq30"));
    for index in [6, 4, 5, 3, 2, 1, 0] {
        let fields: Vec<_> = lines.next().unwrap().split('\t').collect();
        assert_eq!(
            fields[0],
            format!("{}{}", "A".repeat(16), (b'A' + index as u8) as char)
        );
        assert_eq!(fields[1].parse::<u64>().unwrap(), values[index]);
        for (field, value) in [
            (fields[2], values[index] / 3),
            (fields[3], values[index] - values[index] / 3),
        ] {
            assert_eq!(
                field
                    .split(',')
                    .map(|n| n.parse::<u64>().unwrap())
                    .collect::<Vec<_>>(),
                vec![value; 17]
            );
        }
    }
    assert!(lines.next().is_none());
    let metrics: droplet::CountMetrics =
        serde_json::from_slice(&std::fs::read(root.join("metrics.json")).unwrap()).unwrap();
    assert_eq!(metrics.unique_barcodes, 7);
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn count_table_rehash_keeps_all_counts_and_global_order() {
    let root = std::env::temp_dir().join(format!("demux-sharded-counts-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let mut expected = BTreeMap::new();
    let mut actual = CountTable::default();
    for i in 0..20_000u64 {
        let mut key = [b'A'; 17];
        let mut number = i;
        for base in &mut key {
            *base = b"AGT"[(number % 3) as usize];
            number /= 3;
        }
        let value = droplet::BarcodeCount {
            count: i % 1000 + 1,
            below_q20: [i % 300; 17],
            q30: [i % 700; 17],
        };
        actual.merge(key, &value).unwrap();
        actual.merge(key, &value).unwrap();
        expected.insert(
            String::from_utf8(key.to_vec()).unwrap(),
            (value.count * 2, value.below_q20[0] * 2, value.q30[0] * 2),
        );
    }
    assert_eq!(actual.len(), 20_000);
    assert!(actual.parts.iter().all(|part| !part.is_empty()));
    write_compact_counts(
        actual,
        Default::default(),
        &root.join("counts.tsv"),
        &root.join("metrics.json"),
    )
    .unwrap();
    let output = std::fs::read_to_string(root.join("counts.tsv")).unwrap();
    let mut previous: Option<(u64, String)> = None;
    for line in output.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        let (count, low, high) = expected.remove(fields[0]).unwrap();
        assert_eq!(fields[1].parse::<u64>().unwrap(), count);
        for (field, value) in [(fields[2], low), (fields[3], high)] {
            assert_eq!(
                field
                    .split(',')
                    .map(|n| n.parse::<u64>().unwrap())
                    .collect::<Vec<_>>(),
                vec![value; 17]
            );
        }
        if let Some((n, barcode)) = previous {
            assert!(n > count || (n == count && barcode.as_str() < fields[0]));
        }
        previous = Some((count, fields[0].to_owned()));
    }
    assert!(expected.is_empty());
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
#[ignore = "Linux16GiB高基数内存压力验收，普通回归不分配10GiB"]
fn count_table_large_cardinality_under_memory_limit() {
    let mut ballast = vec![0u8; 10usize * 1024 * 1024 * 1024];
    for byte in ballast.iter_mut().step_by(4096) {
        *byte = 1;
    }
    let mut counts = CountTable::default();
    let added = droplet::BarcodeCount {
        count: 1,
        below_q20: [0; 17],
        q30: [1; 17],
    };
    for i in 0..40_000_000u64 {
        let mut key = [b'A'; 17];
        let mut number = i;
        for base in &mut key {
            *base = b"AGT"[(number % 3) as usize];
            number /= 3;
        }
        counts.merge(key, &added).unwrap();
    }
    assert_eq!(counts.len(), 40_000_000);
    assert!(counts
        .iter()
        .all(|(_, value)| value.count == 1 && value.expanded().q30 == [1; 17]));
    std::hint::black_box(&ballast);
    std::hint::black_box(&counts);
    eprintln!(
        "DROPLET_COUNT_MEMORY_STRESS entries={} ballast_bytes={} shards={}",
        counts.len(),
        ballast.len(),
        counts.parts.len()
    );
}

struct Staged {
    parse_seconds: f64,
    compress_seconds: f64,
    input: InputChunk,
    metrics: droplet::CountMetrics,
    counts: AHashMap<[u8; 17], droplet::BarcodeCount>,
    block: Block,
    compressed: Vec<Vec<u8>>,
}
fn stage_batch(id: u64, input: InputChunk, shards: usize) -> Result<Staged> {
    let started = Instant::now();
    let mut raw: Vec<Vec<u8>> = (0..shards).map(|_| Vec::new()).collect();
    let mut counts: AHashMap<[u8; 17], droplet::BarcodeCount> = AHashMap::new();
    let mut metrics = droplet::CountMetrics::default();
    let mut block = Block {
        id,
        pairs: input.count as u64,
        base_ambiguous: 0,
        base_unmatched: 0,
        base_short: 0,
        segments: Vec::new(),
    };
    for i in 0..input.count {
        let seq = &input.r1_seqs[input.r1_offsets[i]..input.r1_offsets[i + 1]];
        let qual = &input.r1_quals[input.r1_offsets[i]..input.r1_offsets[i + 1]];
        let r2 = &input.r2_seqs[input.r2_offsets[i]..input.r2_offsets[i + 1]];
        let q2 = &input.r2_quals[input.r2_offsets[i]..input.r2_offsets[i + 1]];
        anyhow::ensure!(
            seq.len() == qual.len() && r2.len() == q2.len(),
            "sequence/quality mismatch"
        );
        anyhow::ensure!(
            qual.iter().chain(q2).all(|q| (33..=126).contains(q)),
            "invalid quality"
        );
        metrics.input_reads += 1;
        match droplet::extract(seq) {
            Ok(layout) => {
                metrics.structured_reads += 1;
                metrics.barcode_with_c += u64::from(layout.barcode.contains(&b'C'));
                metrics.umi_with_n += u64::from(layout.umi.contains(&b'N'));
                let count = counts.entry(layout.barcode.try_into()?).or_default();
                count.count += 1;
                for (j, q) in qual[layout.barcode_start..layout.barcode_start + 17]
                    .iter()
                    .enumerate()
                {
                    count.below_q20[j] += u64::from(*q < 53);
                    count.q30[j] += u64::from(*q >= 63);
                }
                let id1 = &input.r1_ids[input.r1_id_offsets[i]..input.r1_id_offsets[i + 1]];
                let id2 = &input.r2_ids[input.r2_id_offsets[i]..input.r2_id_offsets[i + 1]];
                let out = &mut raw[shard(layout.barcode, shards)];
                append_u32(out, i)?;
                out.extend_from_slice(layout.barcode);
                out.push(u8::from(r2.len() <= 9));
                let key = droplet::molecule_key(&layout, seq);
                out.push(u8::from(key.is_some()));
                out.extend_from_slice(&key.unwrap_or(0).to_le_bytes());
                let length_at = out.len();
                out.extend_from_slice(&[0; 8]);
                let first = out.len();
                write_fastq_record(
                    out,
                    id1,
                    seq,
                    qual,
                    layout.trim_start,
                    Some(layout.umi),
                    true,
                    false,
                )?;
                let second = out.len();
                write_fastq_record(out, id2, r2, q2, 9, Some(layout.umi), true, false)?;
                let n1 = u32::try_from(second - first)?;
                let n2 = u32::try_from(out.len() - second)?;
                out[length_at..length_at + 4].copy_from_slice(&n1.to_le_bytes());
                out[length_at + 4..length_at + 8].copy_from_slice(&n2.to_le_bytes());
            }
            Err(fate) => {
                match fate {
                    TerminalFate::Ambiguous => metrics.ambiguous += 1,
                    TerminalFate::Short => metrics.short += 1,
                    _ => metrics.unmatched += 1,
                }
                if seq.len() < MIN_R1_LEN || fate == TerminalFate::Short {
                    block.base_short += 1;
                } else if fate == TerminalFate::Ambiguous {
                    block.base_ambiguous += 1;
                } else {
                    block.base_unmatched += 1;
                }
            }
        }
    }
    let parse_seconds = started.elapsed().as_secs_f64();
    let compressed_at = Instant::now();
    let compressed = compress_stage_shards(&raw);
    block.segments = raw
        .iter()
        .map(|x| Segment {
            offset: 0,
            bytes: 0,
            raw_bytes: x.len() as u64,
        })
        .collect();
    let compressed = compressed?;
    Ok(Staged {
        input,
        metrics,
        counts,
        block,
        compressed,
        parse_seconds,
        compress_seconds: compressed_at.elapsed().as_secs_f64(),
    })
}
struct MateBatch {
    ids: Vec<u8>,
    seqs: Vec<u8>,
    quals: Vec<u8>,
    offsets: Vec<usize>,
    id_offsets: Vec<usize>,
    count: usize,
    lane: usize,
}
impl MateBatch {
    fn new() -> Self {
        Self {
            ids: Vec::new(),
            seqs: Vec::new(),
            quals: Vec::new(),
            offsets: vec![0],
            id_offsets: vec![0],
            count: 0,
            lane: 0,
        }
    }
    fn clear(&mut self) {
        self.ids.clear();
        self.seqs.clear();
        self.quals.clear();
        self.offsets.truncate(1);
        self.id_offsets.truncate(1);
        self.count = 0;
    }
    fn swap_into(&mut self, chunk: &mut InputChunk, mate: usize) {
        let (ids, seqs, quals, offsets, id_offsets) = if mate == 1 {
            (
                &mut chunk.r1_ids,
                &mut chunk.r1_seqs,
                &mut chunk.r1_quals,
                &mut chunk.r1_offsets,
                &mut chunk.r1_id_offsets,
            )
        } else {
            (
                &mut chunk.r2_ids,
                &mut chunk.r2_seqs,
                &mut chunk.r2_quals,
                &mut chunk.r2_offsets,
                &mut chunk.r2_id_offsets,
            )
        };
        std::mem::swap(ids, &mut self.ids);
        std::mem::swap(seqs, &mut self.seqs);
        std::mem::swap(quals, &mut self.quals);
        std::mem::swap(offsets, &mut self.offsets);
        std::mem::swap(id_offsets, &mut self.id_offsets);
        self.clear();
    }
}
struct DecoderChild(Option<std::process::Child>);
impl Drop for DecoderChild {
    fn drop(&mut self) {
        if let Some(child) = &mut self.0 {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}
fn read_mate(
    paths: Vec<PathBuf>,
    mate: usize,
    batch: usize,
    decoder_python: PathBuf,
    decoder_threads: usize,
    tx: Sender<MateBatch>,
    recycle: Receiver<MateBatch>,
) -> Result<u64> {
    let started = Instant::now();
    let mut send_seconds = 0.0;
    let mut parse_seconds = 0.0;
    let mut total = 0u64;
    for (lane, path) in paths.into_iter().enumerate() {
        ensure_fastq_input(&path, if mate == 1 { "R1" } else { "R2" })?;
        let mut child = DecoderChild(None);
        let decoder_threads = decoder_threads.to_string();
        let mut process = std::process::Command::new(&decoder_python)
            .args([
                "-c",
                "import rapidgzip; raise SystemExit(rapidgzip.cli())",
                "-P",
                &decoder_threads,
                "--io-read-method",
                "sequential",
                "--verify",
                "-d",
                "-c",
            ])
            .arg(&path)
            .stdout(std::process::Stdio::piped())
            .spawn()?;
        let output = process.stdout.take().context("missing rapidgzip stdout")?;
        child.0 = Some(process);
        let mut reader = needletail::parse_fastx_reader(output)?;
        let mut ended = false;
        while !ended {
            let mut chunk = recycle.try_recv().unwrap_or_else(|_| MateBatch::new());
            chunk.clear();
            chunk.lane = lane;
            let at = Instant::now();
            while chunk.count < batch {
                match reader.next() {
                    Some(record) => {
                        let record = record?;
                        let seq = record.seq();
                        let qual = record.qual().context("missing FASTQ quality")?;
                        anyhow::ensure!(seq.len() == qual.len(), "sequence/quality mismatch");
                        chunk.ids.extend_from_slice(record.id());
                        chunk.seqs.extend_from_slice(seq.as_ref());
                        chunk.quals.extend_from_slice(qual);
                        chunk.offsets.push(chunk.seqs.len());
                        chunk.id_offsets.push(chunk.ids.len());
                        chunk.count += 1;
                        total += 1;
                    }
                    None => {
                        ended = true;
                        break;
                    }
                }
            }
            parse_seconds += at.elapsed().as_secs_f64();
            if chunk.count > 0 {
                let at = Instant::now();
                tx.send(chunk)?;
                send_seconds += at.elapsed().as_secs_f64();
            }
        }
        drop(reader);

        if let Some(mut process) = child.0.take() {
            anyhow::ensure!(
                process.wait()?.success(),
                "rapidgzip failed for {}",
                path.display()
            );
        }
    }
    eprintln!(
        "DROPLET_READER_PROFILE {}",
        serde_json::json!({"mode":"parallel_mate","mate":mate,"pairs":total,"wall_seconds":started.elapsed().as_secs_f64(),"parse_seconds":parse_seconds,"send_wait_seconds":send_seconds})
    );
    Ok(total)
}
fn produce_parallel(
    args: &StageArgs,
    recycle: Receiver<InputChunk>,
    tx: Sender<(u64, InputChunk)>,
) -> Result<u64> {
    let (tx1, rx1) = bounded(2);
    let (tx2, rx2) = bounded(2);
    let (pool1, free1) = crossbeam_channel::unbounded();
    let (pool2, free2) = crossbeam_channel::unbounded();
    let r1 = args.r1.clone();
    let r2 = args.r2.clone();
    let batch = args.batch_pairs;
    let rapid1 = args.decoder_python.clone();
    let rapid2 = args.decoder_python.clone();
    let decoder_threads = (args.threads / 4).clamp(1, 4);
    let first =
        std::thread::spawn(move || read_mate(r1, 1, batch, rapid1, decoder_threads, tx1, free1));
    let second =
        std::thread::spawn(move || read_mate(r2, 2, batch, rapid2, decoder_threads, tx2, free2));
    let mut total = 0u64;
    let mut block = 0u64;
    let started = Instant::now();
    let mut receive = 0.0;
    let mut downstream = 0.0;
    let result = (|| -> Result<()> {
        loop {
            let at = Instant::now();
            let a = rx1.recv();
            let b = rx2.recv();
            receive += at.elapsed().as_secs_f64();
            match (a, b) {
                (Ok(mut a), Ok(mut b)) => {
                    anyhow::ensure!(
                        a.lane == b.lane && a.count == b.count,
                        "unpaired mate lane/batch"
                    );
                    for i in 0..a.count {
                        validate_paired_read_ids(
                            &a.ids[a.id_offsets[i]..a.id_offsets[i + 1]],
                            &b.ids[b.id_offsets[i]..b.id_offsets[i + 1]],
                        )?;
                    }
                    let at = Instant::now();
                    let mut chunk = recycle.recv()?;
                    chunk.count = a.count;
                    total += a.count as u64;
                    a.swap_into(&mut chunk, 1);
                    b.swap_into(&mut chunk, 2);
                    let _ = pool1.send(a);
                    let _ = pool2.send(b);
                    tx.send((block, chunk))?;
                    block += 1;
                    downstream += at.elapsed().as_secs_f64();
                }
                (Err(_), Err(_)) => break,
                _ => anyhow::bail!("unpaired mate EOF"),
            }
        }
        Ok(())
    })();
    drop(rx1);
    drop(rx2);
    drop(pool1);
    drop(pool2);
    let left = first
        .join()
        .map_err(|_| anyhow::anyhow!("R1 reader panic"))?;
    let right = second
        .join()
        .map_err(|_| anyhow::anyhow!("R2 reader panic"))?;
    result?;
    anyhow::ensure!(left? == total && right? == total, "mate total mismatch");
    eprintln!(
        "DROPLET_PAIR_PROFILE {}",
        serde_json::json!({"pairs":total,"blocks":block,"wall_seconds":started.elapsed().as_secs_f64(),"receive_wait_seconds":receive,"downstream_wait_seconds":downstream})
    );
    Ok(total)
}

pub(crate) fn run() -> Result<()> {
    let mut argv: Vec<_> = std::env::args_os().collect();
    argv.remove(1);
    let args = StageArgs::parse_from(argv);
    anyhow::ensure!(
        args.threads > 0 && args.shards > 0 && args.batch_pairs > 0,
        "invalid stage budgets"
    );
    anyhow::ensure!(
        !args.r1.is_empty() && args.r1.len() == args.r2.len(),
        "paired input list mismatch"
    );
    let inputs: Vec<_> = args.r1.iter().map(|p| ("--r1", p.as_path()))
        .chain(args.r2.iter().map(|p| ("--r2", p.as_path())))
        .chain(std::iter::once(("--decoder-python", args.decoder_python.as_path())))
        .collect();
    crate::cli::validate_io_paths(&inputs, &[("--spool", args.spool.as_path())],
        &[("--counts-output", args.counts_output.as_path()),
          ("--count-metrics", args.count_metrics.as_path())])?;
    let compute_threads = args.threads.saturating_sub(13).max(1);
    rayon::ThreadPoolBuilder::new()
        .num_threads(compute_threads)
        .build_global()?;
    std::fs::create_dir(&args.spool)?;
    let capacity = compute_threads.saturating_add(1).max(2);
    let (tx_recycle, rx_recycle) = bounded(capacity);
    let (tx_work, rx_work) = bounded(capacity);
    let (tx_done, rx_done) = bounded::<(u64, Staged)>(capacity);
    for _ in 0..capacity {
        tx_recycle.send(InputChunk::new(args.batch_pairs))?;
    }
    let shards = args.shards;
    let reader_args = StageArgs::parse_from({
        let mut a: Vec<_> = std::env::args_os().collect();
        a.remove(1);
        a
    });
    let reader = std::thread::spawn(move || produce_parallel(&reader_args, rx_recycle, tx_work));
    let worker = std::thread::spawn(move || -> Result<()> {
        rx_work
            .into_iter()
            .par_bridge()
            .try_for_each_with(tx_done, |tx, (id, input)| {
                tx.send((id, stage_batch(id, input, shards)?))?;
                Ok::<(), anyhow::Error>(())
            })?;
        Ok(())
    });
    let (tx_spool, rx_spool) = bounded::<(Block, Vec<Vec<u8>>, InputChunk)>(2);
    let spool_path = args.spool.clone();
    let spool_recycle = tx_recycle.clone();
    let spool_writer = std::thread::spawn(move || -> Result<u64> {
        let writeback = writeback::WritebackGuard::start(DemuxMode::DnaOnlyDroplet, shards)?;
        let mut writers: Vec<_> = (0..shards)
            .map(|i| {
                let file = File::create(spool_path.join(format!("{i:04}.pack")))?;
                Ok(BufWriter::with_capacity(
                    1024 * 1024,
                    writeback::spool_output(file, writeback.as_ref())?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut index = BufWriter::new(File::create(spool_path.join("index.jsonl"))?);
        let mut offsets = vec![0u64; shards];
        let mut next = 0u64;
        let mut write_seconds = 0.0;
        for (mut block, data, mut input) in rx_spool {
            anyhow::ensure!(block.id == next, "spool block order mismatch");
            next += 1;
            let at = Instant::now();
            for (i, data) in data.into_iter().enumerate() {
                block.segments[i].offset = offsets[i];
                block.segments[i].bytes = data.len() as u64;
                writers[i].write_all(&data)?;
                offsets[i] += data.len() as u64;
            }
            serde_json::to_writer(&mut index, &block)?;
            index.write_all(b"\n")?;
            write_seconds += at.elapsed().as_secs_f64();
            input.clear();
            let _ = spool_recycle.send(input);
        }
        let at = Instant::now();
        for writer in &mut writers {
            writer.flush()?;
            writer.get_ref().sync_spool()?;
        }
        if let Some(guard) = writeback {
            guard.finish()?;
        }
        index.flush()?;
        index.get_ref().sync_all()?;
        eprintln!(
            "DROPLET_SPOOL_WRITER_PROFILE {}",
            serde_json::json!({"blocks":next,"write_seconds":write_seconds,"flush_sync_seconds":at.elapsed().as_secs_f64(),"spool_bytes":offsets.iter().sum::<u64>()})
        );
        Ok(next)
    });
    let mut fixed_counts = CountTable::default();
    let mut metrics = droplet::CountMetrics::default();
    let mut pending = BTreeMap::new();
    let mut next = 0u64;
    let start = Instant::now();
    let mut worker_parse = 0.0;
    let mut worker_compress = 0.0;
    let mut merge_seconds = 0.0;
    let mut write_seconds = 0.0;
    let mut receive_seconds = 0.0;
    let consumed = (|| -> Result<()> {
        while let Ok((id, item)) = {
            let at = Instant::now();
            let result = rx_done.recv();
            receive_seconds += at.elapsed().as_secs_f64();
            result
        } {
            for item in insert_completed_chunk(&mut pending, &mut next, id, item)? {
                worker_parse += item.parse_seconds;
                worker_compress += item.compress_seconds;
                let at = Instant::now();
                for (bc, stat) in item.counts {
                    fixed_counts.merge(bc, &stat)?;
                }
                merge_seconds += at.elapsed().as_secs_f64();
                metrics.input_reads += item.metrics.input_reads;
                metrics.structured_reads += item.metrics.structured_reads;
                metrics.ambiguous += item.metrics.ambiguous;
                metrics.unmatched += item.metrics.unmatched;
                metrics.short += item.metrics.short;
                metrics.barcode_with_c += item.metrics.barcode_with_c;
                metrics.umi_with_n += item.metrics.umi_with_n;
                let at = Instant::now();
                tx_spool.send((item.block, item.compressed, item.input))?;
                write_seconds += at.elapsed().as_secs_f64();
                if next % 25 == 0 {
                    log::info!(
                        "raw stage pairs={} wall={:.3}",
                        metrics.input_reads,
                        start.elapsed().as_secs_f64()
                    );
                }
            }
        }
        Ok(())
    })();
    drop(rx_done);
    drop(tx_spool);
    drop(tx_recycle);
    let written = spool_writer
        .join()
        .map_err(|_| anyhow::anyhow!("spool writer panic"));
    let read_result = reader.join().map_err(|_| anyhow::anyhow!("reader panic"));
    let work_result = worker.join().map_err(|_| anyhow::anyhow!("worker panic"));
    anyhow::ensure!(written?? == next, "spool block count mismatch");
    consumed?;
    anyhow::ensure!(pending.is_empty(), "missing staged block");
    let total = read_result??;
    work_result??;
    anyhow::ensure!(total == metrics.input_reads, "stage count mismatch");
    let at = Instant::now();
    write_compact_counts(
        fixed_counts,
        metrics,
        &args.counts_output,
        &args.count_metrics,
    )?;
    eprintln!(
        "DROPLET_STAGE_PROFILE {}",
        serde_json::json!({"compute_threads":compute_threads,"inflight_chunks":capacity,"inline_count_bytes":std::mem::size_of::<CompactCount>(),"count_table_shards":64,"batch_pairs":args.batch_pairs,"rank_key_bytes":std::mem::size_of::<([u8;17],u64)>(),"pairs":total,"blocks":next,"wall_seconds":start.elapsed().as_secs_f64(),"worker_parse_seconds_sum":worker_parse,"worker_compress_seconds_sum":worker_compress,"main_count_merge_seconds":merge_seconds,"main_spool_queue_wait_seconds":write_seconds,"main_receive_wait_seconds":receive_seconds,"final_counts_seconds":at.elapsed().as_secs_f64()})
    );
    let manifest = Manifest {
        version: 1,
        shards,
        blocks: next,
        pairs: total,
        input_fastq_pairs: args.r1.len(),
    };
    let temp = args.spool.join("manifest.json.tmp");
    std::fs::write(&temp, serde_json::to_vec_pretty(&manifest)?)?;
    File::open(&temp)?.sync_all()?;
    std::fs::rename(temp, args.spool.join("manifest.json"))?;
    log::info!(
        "raw stage complete pairs={total} blocks={next} wall={:.3}",
        start.elapsed().as_secs_f64()
    );
    Ok(())
}

//! Stage写出的共享格式与严格解码；replay和分子evidence消费同一校验。
use crate::*;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};

#[derive(Serialize, serde::Deserialize, Clone)]
pub(super) struct Segment {
    pub(super) offset: u64,
    pub(super) bytes: u64,
    pub(super) raw_bytes: u64,
}
#[derive(Serialize, serde::Deserialize, Clone)]
pub(super) struct Block {
    pub(super) id: u64,
    pub(super) pairs: u64,
    pub(super) base_ambiguous: u64,
    pub(super) base_unmatched: u64,
    pub(super) base_short: u64,
    pub(super) segments: Vec<Segment>,
}
#[derive(Serialize, serde::Deserialize)]
pub(super) struct Manifest {
    pub(super) version: u32,
    pub(super) shards: usize,
    pub(super) blocks: u64,
    pub(super) pairs: u64,
    pub(super) input_fastq_pairs: usize,
}
pub(super) fn shard(barcode: &[u8], n: usize) -> usize {
    let mut h = 0xcbf29ce484222325u64;
    for base in barcode {
        h = (h ^ u64::from(*base)).wrapping_mul(0x100000001b3);
    }
    (h as usize) % n
}
pub(super) fn append_u32(out: &mut Vec<u8>, n: usize) -> Result<()> {
    out.extend_from_slice(&u32::try_from(n)?.to_le_bytes());
    Ok(())
}
pub(super) struct RawRecord<'a> {
    pub(super) ordinal: u32,
    pub(super) barcode: &'a [u8],
    pub(super) short_r2: bool,
    pub(super) key: Option<u128>,
    pub(super) r1: &'a [u8],
    pub(super) r2: &'a [u8],
}

// CRC 只证明压缩字节完整；重放和分子证据还须核对内部 FASTQ 结构与身份。
pub(super) fn spool_fastq(data: &[u8]) -> Result<(&[u8], &[u8])> {
    let mut lines = data.split(|b| *b == b'\n');
    let id = lines.next().context("missing spool FASTQ header")?;
    let seq = lines.next().context("missing spool FASTQ sequence")?;
    let plus = lines.next().context("missing spool FASTQ separator")?;
    let qual = lines.next().context("missing spool FASTQ quality")?;
    anyhow::ensure!(
        id.starts_with(b"@")
            && id.len() > 1
            && !id.iter().any(u8::is_ascii_whitespace)
            && plus == b"+"
            && !seq.is_empty()
            && lines.next() == Some(b"".as_slice())
            && lines.next().is_none(),
        "invalid spool FASTQ structure"
    );
    validate_fastq_quality(seq, qual)?;
    Ok((&id[1..], seq))
}

pub(super) fn walk(
    data: &[u8],
    mut consume: impl FnMut(RawRecord<'_>) -> Result<()>,
) -> Result<()> {
    let mut at = 0usize;
    while at < data.len() {
        anyhow::ensure!(data.len() - at >= 47, "truncated record header");
        let ordinal = u32::from_le_bytes(data[at..at + 4].try_into()?);
        let barcode = &data[at + 4..at + 21];
        anyhow::ensure!(
            barcode.iter().all(|b| b"ACGT".contains(b)),
            "invalid spool barcode"
        );
        anyhow::ensure!(
            data[at + 21] <= 1 && data[at + 22] <= 1,
            "invalid spool boolean flag"
        );
        let short_r2 = data[at + 21] != 0;
        let encoded = u128::from_le_bytes(data[at + 23..at + 39].try_into()?);
        let key = (data[at + 22] != 0).then_some(encoded);
        anyhow::ensure!(
            key.is_some() || encoded == 0,
            "invalid absent spool molecule key"
        );
        let n1 = u32::from_le_bytes(data[at + 39..at + 43].try_into()?) as usize;
        let n2 = u32::from_le_bytes(data[at + 43..at + 47].try_into()?) as usize;
        at += 47;
        let end = at
            .checked_add(n1)
            .and_then(|x| x.checked_add(n2))
            .context("record length overflow")?;
        anyhow::ensure!(end <= data.len(), "truncated record body");
        let r1 = &data[at..at + n1];
        let r2 = &data[at + n1..end];
        let (id1, seq1) = spool_fastq(r1)?;
        let (id2, _) = spool_fastq(r2)?;
        anyhow::ensure!(id1 == id2, "spool mate QNAME mismatch");
        let umi_start = id1
            .iter()
            .rposition(|b| *b == b':')
            .context("missing spool UMI")?;
        let umi = &id1[umi_start + 1..];
        anyhow::ensure!(
            umi_start > 0 && umi.len() == 12 && umi.iter().all(|b| b"ACGTN".contains(b)),
            "invalid spool UMI"
        );
        let expected = seq1
            .get(..32)
            .and_then(|insert| droplet::encode_molecule(umi, insert));
        anyhow::ensure!(key == expected, "spool molecule key disagrees with FASTQ");
        consume(RawRecord {
            ordinal,
            barcode,
            short_r2,
            key,
            r1,
            r2,
        })?;
        at = end;
    }
    Ok(())
}
pub(super) fn blocks(spool: &Path) -> Result<(Manifest, Vec<Block>)> {
    let manifest: Manifest = serde_json::from_slice(&std::fs::read(spool.join("manifest.json"))?)?;
    anyhow::ensure!(
        manifest.version == 1 && manifest.shards > 0,
        "unsupported spool"
    );
    let mut out = Vec::new();
    for line in BufReader::new(File::open(spool.join("index.jsonl"))?).lines() {
        let block: Block = serde_json::from_str(&line?)?;
        anyhow::ensure!(
            block.id == out.len() as u64 && block.segments.len() == manifest.shards,
            "invalid block sequence"
        );
        out.push(block);
    }
    anyhow::ensure!(
        out.len() as u64 == manifest.blocks
            && out.iter().map(|b| b.pairs).sum::<u64>() == manifest.pairs,
        "invalid spool totals"
    );
    for shard in 0..manifest.shards {
        let mut offset = 0;
        for block in &out {
            anyhow::ensure!(
                block.segments[shard].offset == offset,
                "noncontiguous shard"
            );
            offset += block.segments[shard].bytes;
        }
        anyhow::ensure!(
            std::fs::metadata(spool.join(format!("{shard:04}.pack")))?.len() == offset,
            "shard size mismatch"
        );
    }
    Ok((manifest, out))
}
pub(super) fn decode(spool: &Path, shard: usize, segment: &Segment) -> Result<Vec<u8>> {
    if segment.bytes == 0 {
        anyhow::ensure!(segment.raw_bytes == 0, "empty compressed segment");
        return Ok(Vec::new());
    }
    anyhow::ensure!(
        segment.raw_bytes <= 1024 * 1024 * 1024 && segment.bytes <= 1024 * 1024 * 1024,
        "oversized segment"
    );
    let mut f = File::open(spool.join(format!("{shard:04}.pack")))?;
    f.seek(SeekFrom::Start(segment.offset))?;
    let mut gz = vec![0u8; segment.bytes as usize];
    f.read_exact(&mut gz)?;
    let mut raw = vec![0u8; segment.raw_bytes as usize];
    let n = libdeflate::Decompressor::new().gzip_decompress(&gz, &mut raw)?;
    anyhow::ensure!(n == raw.len(), "decompressed size mismatch");
    Ok(raw)
}

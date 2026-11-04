//! Stage打包与最终配对输出共用的压缩实现；必须进入calling源码身份。
use anyhow::{Context, Result};
use libdeflate::{CompressionLvl, Compressor};
use rayon::prelude::*;

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

pub(super) fn compress_stage_shards(raw: &[Vec<u8>]) -> Result<Vec<Vec<u8>>> {
    raw.par_iter()
        .map(|plain| {
            let mut out = Vec::new();
            if !plain.is_empty() {
                let mut compressor = Compressor::new(CompressionLvl::new(1).unwrap());
                compress_plain_buffer(&mut compressor, plain, &mut out, "raw barcode shard")?;
            }
            Ok(out)
        })
        .collect()
}

//! 直接扫描完整spool核验独立分子证据，不以父条码诊断代替观测。
use super::spool::*;
use crate::*;

pub(crate) fn run(args: &Args, spool: &Path) -> Result<()> {
    #[derive(serde::Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Pair {
        barcode: String,
        parent: String,
    }
    let pairs: Vec<Pair> =
        serde_json::from_slice(&std::fs::read(args.evidence_pairs.as_ref().unwrap())?)?;
    anyhow::ensure!(!pairs.is_empty(), "empty evidence request");
    let (manifest, blocks) = blocks(spool)?;
    let mut children: AHashMap<Vec<u8>, droplet::Molecules> = AHashMap::new();
    for pair in &pairs {
        anyhow::ensure!(
            [&pair.barcode, &pair.parent]
                .into_iter()
                .all(|bc| bc.len() == 17 && bc.bytes().all(|b| b"AGT".contains(&b)))
                && pair.barcode != pair.parent,
            "Invalid molecule evidence pair"
        );
        anyhow::ensure!(
            children
                .insert(
                    pair.barcode.as_bytes().to_vec(),
                    droplet::Molecules::default()
                )
                .is_none(),
            "duplicate child"
        );
    }
    let scan = |selected: &AHashSet<Vec<u8>>, consume: &mut dyn FnMut(&[u8], u128)| -> Result<()> {
        let wanted: AHashSet<_> = selected
            .iter()
            .map(|bc| shard(bc, manifest.shards))
            .collect();
        for block in &blocks {
            for shard in &wanted {
                let raw = decode(spool, *shard, &block.segments[*shard])?;
                walk(&raw, |record| {
                    if selected.contains(record.barcode) {
                        if let Some(key) = record.key {
                            consume(record.barcode, key);
                        }
                    }
                    Ok(())
                })?;
            }
        }
        Ok(())
    };
    scan(&children.keys().cloned().collect(), &mut |bc, key| {
        children.get_mut(bc).unwrap().add(key)
    })?;
    let mut targets: AHashMap<Vec<u8>, AHashMap<u128, Vec<usize>>> = AHashMap::new();
    let mut parent_reads: AHashMap<Vec<u8>, u64> = AHashMap::new();
    let mut shared: Vec<AHashSet<u128>> = pairs.iter().map(|_| AHashSet::new()).collect();
    for (i, pair) in pairs.iter().enumerate() {
        let lookup = targets.entry(pair.parent.as_bytes().to_vec()).or_default();
        parent_reads
            .entry(pair.parent.as_bytes().to_vec())
            .or_default();
        for (_, key) in &children[pair.barcode.as_bytes()].keys {
            lookup.entry(*key).or_default().push(i);
        }
    }
    scan(&targets.keys().cloned().collect(), &mut |bc, key| {
        *parent_reads.get_mut(bc).unwrap() += 1;
        if let Some(ids) = targets[bc].get(&key) {
            for id in ids {
                shared[*id].insert(key);
            }
        }
    })?;
    let evidence:BTreeMap<_,_>=pairs.iter().enumerate().map(|(i,pair)| {
            let child=&children[pair.barcode.as_bytes()]; let mut keys:Vec<_>=shared[i].iter().copied().collect(); keys.sort_unstable();
            (pair.barcode.clone(),serde_json::json!({"parent":pair.parent,"child_eligible_reads":child.reads,
                "child_sampled_signatures":child.keys.len(),"child_saturated":child.saturated,
                "parent_eligible_reads":parent_reads[pair.parent.as_bytes()],
                "shared_signatures":keys.into_iter().map(|key|format!("{key:022x}")).collect::<Vec<_>>()}))
        }).collect();
    std::fs::write(
        args.evidence_output.as_ref().unwrap(),
        serde_json::to_vec(&serde_json::json!({
            "signature":"exact_umi12_insert32","child_cap":4096,"parent_scan":"complete","pairs":evidence}))?,
    )?;
    Ok(())
}

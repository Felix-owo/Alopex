# Alopex

English · [简体中文](README.zh-CN.md)

Alopex turns paired-end single-cell DNA methylation FASTQs into per-cell CpG counts, conversion and alignment QC, and a sealed result inventory. For a configured project, run `./run_pipeline.sh --dry-run` from the project directory, check the selected libraries and backend, then run `./run_pipeline.sh`.

[Choose a protocol](#choose-a-protocol) · [Install and check](#install-and-check) · [Configure a project](#configure-a-project) · [Run and recover](#run-and-recover) · [Read results](#read-results)

## Choose a protocol

Set `analysis.protocol` once per project and use the same spelling in its manifest. Alopex v19 supports both DD chemistries.

| Protocol | DNA processing | Matching RNA |
| --- | --- | --- |
| `DD-MET3` | R1 contains CB17, 7F, 17L, ME, gap9 and the insert. The first 12 bases of the trimmed insert form the fragment tag and remain in the alignment input. | Lemmus DD-MET3, 3′ DDV2. Match the official RNA/DNA barcode pair within the linked library. |
| `DD-MET5` | R1 contains CB17, the library UMI12, TSO, 17L, ME, gap9 and the insert. Carry the original UMI through alignment and cytosine consensus. | Lemmus DD-MET5, 5′ RNA. Match the same CB17 within the linked library. |
| `cabernet` | Plate barcode map, Bismark by default or optional BISCUIT, paired-end duplicate handling and high-CpH filtering. | Explicit library and plate identities. |
| `taps` | Plate barcode map, BWA-MEM and Rastair. Publishes CpG counts, SNP calls and non-CpG QC. | Optional external RNA library link. |

Both DD protocols independently call DNA cells using raw structured-read counts, a histogram valley and 829,440 DNA design barcodes. DD-MET3 records suspect barcode quality as a diagnostic. DD-MET5 can remove a suspect barcode with independent evidence from library UMIs and insert sequences. Both trim nine bases from R2, use non-directional Bismark and count consensus observations by cell, strand and reference cytosine. Bases with unresolved ties contribute no observation. The BAM retains PCR copies. The DNA design is shared across both chemistries. [Library structures](https://github.com/seekgene/SeekSoulMethyl/tree/v2.2.1#data-structure) · [Official fragment-tag implementation](https://github.com/seekgene/SeekSoulMethyl/blob/v2.2.1/nf/bin/barcode_cs_multi.py#L1219-L1233)

## Install and check

1. Install Conda, Git, curl, gzip and a C compiler. Keep the Pipeline checkout separate from projects. Doctor builds environments, Rust and reference indexes for the target platform. Download `Alopex-v19.2.tar.gz`, `SHA256SUMS` and `SOURCE_SHA256SUMS`, then verify the archive before extracting it. [v19.2 release](https://github.com/Felix-owo/Alopex/releases/tag/v19.2)
2. Uncompress the six source files into `resources/reference_source/` using the names below. Doctor builds the analysis FASTAs, adds lambda and pUC19, builds indexes and records reference identities. The input sources and masked regions are specified in the [reference guide](README.zh-CN.md#2-准备标准-reference).
3. Run the checks below from the Pipeline root. On macOS use the local terminal. On HPC verify and transfer the sealed references first, prepare dependencies on the login node, then run the full Doctor in the terminal of the existing Jupyter compute allocation. Transfer the index families and their producer receipts together. The target platform needs its own managed environments and binaries.
4. Continue after `FINAL STATUS: READY`. Doctor reports its log path under `logs/doctor/`. The standard small scientific runs exercise human Cabernet/Bismark, mouse Cabernet/BISCUIT and TAPS/Rastair.

| Filename | Source |
| --- | --- |
| `hg38.analysis.fa` | [GIABv3 GRCh38 masked no-ALT assembly with decoys](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz) |
| `mm10.fa` | [UCSC GRCm38/mm10](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz) |
| `lambda.fa`, `pUC19.fa` | [Lambda J02459.1](https://www.ncbi.nlm.nih.gov/nuccore/J02459.1) · [pUC19 L09137.2](https://www.ncbi.nlm.nih.gov/nuccore/L09137.2) |
| `hg38.genes.gtf`, `mm10.genes.gtf` | [GENCODE human v48](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz) · [GENCODE mouse M25](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz) |

```bash
bash core/doctor.sh conda
bash core/doctor.sh
```

The managed environment pins Rust Bismark 3.4.0 with SAMtools/HTSlib 1.25. After updating the source, run `bash core/doctor.sh conda` and verify the active binary with `conda/current/bismark/bin/bismark --version`.

To check an installed component, run `bash core/doctor.sh conda --check`, `bash core/doctor.sh demux --check` or `bash core/doctor.sh reference --check`. To rebuild the native demultiplexer, run `bash core/doctor.sh demux`.

Doctor checks FASTA index geometry and the bases of both control accessions. For TSS, the CHR and PRI comprehensive GTFs in the linked human and mouse releases have identical annotations on the primary nuclear chromosomes. Use a single GTF per species. Doctor validates a private replacement before installing a repaired source index. A failed rebuild preserves the previous index.

For a transferred reference, first verify the FASTA and its index receipts in the candidate directory. An explicit import checks the actual file hashes even when the recorded file metadata match. This example checks the two components used by DD-MET3 and DD-MET5. Verify the remaining installed index families and downstream targets in the same way, then run Doctor and a small alignment before replacing a working installation.

```bash
ref=/path/to/candidate/hg38.primary.lambda.puc19.fa
for kind in fasta bismark-combined; do
  PYTHONPATH=core conda/current/control/bin/python -m alopex reference-binding import \
    --reference "$ref" --kind "$kind"
done
```

## Configure a project

```bash
mkdir -p /path/to/Patient001
cd /path/to/Patient001
/path/to/Alopex/core/run_pipeline.sh --init-project
```

Edit the generated `00_config/config.yaml`. For a human DD-MET3 library, set `species` to `hg38` and edit these keys in the full schema-3 template.

```yaml
analysis:
  protocol: DD-MET3
  cpg_coordinate_base: 1
DD-MET3: {}
DD-MET5: {}
```

For DD-MET5 select `DD-MET5`. DD backend, spacer length and library direction are fixed. Keep tunable defaults under `common` or override supported fields in the selected module. For a legacy DD-MET5 configuration, replace `droplet` in `analysis.protocol`, the module name and manifest protocol before running. Existing sealed analyses retain their original software identity. Plate protocols require `00_config/Barcode_Map.csv` with unique `DNA_Barcode,PlateID,Cell_Order` identities. Keep leading zeros in plate and order fields. Use CSV quoting for order labels containing tabs, quotes or line breaks, and double each literal quote within a quoted field.

1. Place paired `.fastq.gz` or `.fq.gz` files under `01_raw/`. Examples include `Lib_R1.fastq.gz`, `Lib_L001_R1_001.fastq.gz` and `Lib.R1.raw.fastq.gz` with their matching R2 files. Put same-named pairs from separate batches of one library in separate subdirectories. Both mates must have the same lane or chunk identity.
2. Run `./run_pipeline.sh --refresh-manifest` to register each library once and report the number of paired chunks. The command returns as soon as the manifest is ready. Review `00_config/sample_manifest.tsv` and keep the intended libraries. Its fields are `sample_id,dna_raw_sample,rna_sample,species,protocol,notes`. Supply the exact Lemmus raw basename in `rna_sample` to link an RNA library. An empty field leaves the DNA library unlinked.
3. Analysis processes all original chunks of each selected library together. At launch, explicit lanes and chunks are ordered by number and same-named batches by directory path. Keep those files stable until delivery. Each physical file can occupy one place in the selected libraries. When a read name ends in `/1` or `/2`, the suffix must match its input side. SeekGene `ID/1:CB17` and `ID/2:CB17` must also agree on the barcode.

## Run and recover

```bash
./run_pipeline.sh --dry-run
./run_pipeline.sh
```

The launcher finds the current project, seals its configuration and inputs and selects the executor. macOS runs locally. Linux/HPC submits rules through Slurm. A dry run checks planning. During a real run, inspect `04_logs/controller/`, `04_logs/rules/` and `04_logs/slurm/`. Project-local `05_tmp/` must have room for the requested scratch space and its margin.

Keep the reference and the selected aligner's indexes stable through publication. A change to a raw FASTQ's ctime, including a permission change, can trigger reprocessing for that library. If the launcher reports changed inputs, verify the original FASTQs and active reference before restarting.

On failure, check the rule and attempt log, keep `02_work/` and `.snakemake/`, fix the cause and rerun. Remove a stale lock only after its controller has stopped. For a failed publication, keep `02_work/results_stage/` and `03_results/` on the same filesystem. Republish from an intact stage. If both paths share a damaged inode, restore verified bytes or regenerate the file with its producing rule. The manifest commits last. A later run retries pending demux cleanup.

If the DNA barcode design fails its hash check, restore `resources/droplet_ME5_U3CB_methylation.txt.gz` from the verified checkout, rebuild with `bash core/doctor.sh demux` and rerun the project. If calling reports no two separated peaks, check the kit and the input read lengths against the read structure before resuming. Read the failed attempt in `04_logs/rules/droplet_calling/<library>.log`.

## Read results

| Path | Contents |
| --- | --- |
| `03_results/run_manifest.json` | Complete delivery marker and file inventory. Verify the snapshot identity and output sizes/hashes before consuming results. |
| `03_results/CpG/<Sample_ID>.cpg.tsv.zst` | Four columns, `chrom, pos, methyl, unmethyl`. The last two are integer counts, and the sealed configuration sets the coordinate base. |
| `03_results/QC_Results/` | Sealed cell manifest, read and conversion metrics, and `multiqc_report.html`. DD cell IDs are `<sample_id>_<DNA barcode>` with an empty plate ID. |
| `03_results/SNP/` | TAPS VCF and index, or optional BISCUIT SNP BED and index. |
| `06_downstream/` | DNA QC, SingleCpG, the DD 1 Mb integer counts for doublet analysis and plots. |

The sealed `run_manifest.json` records the SHA256 of the DNA design. DD QC records `tag_source` as `trimmed_r1_insert_prefix12` for DD-MET3 or `library_umi12` for DD-MET5. Pair duplicate rate is undefined for both. Methylated and unmethylated counts are observations at a particular cytosine on one strand. Keep observed zeros distinct from missing coverage. [Downstream guide](downstream/README.md)

The DD M-bias plot reports methylation by sequencing cycle for physical FASTQ R1 and R2, before tag consensus. To refresh an existing project after updating Alopex, complete Doctor and run `./run_pipeline.sh` from the project directory.

To generate DNA QC, SingleCpG and the 1 Mb doublet input from a sealed DD project, use the managed Notebook interpreter. For plots, select the Alopex QC kernel and run the Notebook from its first cell. Restart the kernel after a checkout update. [Kernel registration](downstream/README.md#jupyter--vs-code-kernel) · [DNA QC Notebook](downstream/DNA_QC_Visualization.ipynb)

```bash
PIPELINE=/path/to/Alopex
PATH="$PIPELINE/conda/current/control/bin:$PATH" PYTHONPATH="$PIPELINE/core" \
  "$PIPELINE/conda/current/notebook/bin/python" -m alopex.downstream.qc \
  --project-dir /path/to/Patient001
```

For joint QC, run the matching Lemmus chemistry on its own reads, then integrate its finalized RNA with the sealed DNA delivery and its current QC. DD-MET3 uses 248,832 official barcode pairs, including 82,944 RNA barcodes outside the RNA quantification whitelist. The pairing table is used only for linking cells called on each side. Both native IDs and the single-modality records are retained.

Software terms are in [LICENSE](LICENSE).

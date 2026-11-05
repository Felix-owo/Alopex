"""Alopex内部命令的参数解析、分派与退出处理。"""
from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

from .config import (
    BISCUIT_INDEX_SUFFIXES,
    BISCUIT_TRANSIENT_SUFFIXES,
    MANAGED_SPECIES,
    SPIKE_CONTIGS,
    backend_active_index_paths,
    envs_section_sha256,
    envs_spec_logical_name,
    export_config_scalars,
    load_project_config,
    managed_alignment_contigs,
    validate_downstream_reference,
    validate_reference_component,
)
from .delivery import (
    create_run_snapshot,
    project_status,
    write_delivery_ready,
    write_final_sample_manifest,
)
from .inputs import (
    _resolve_config,
    call_droplet_cells,
    collect_rna_fastq,
    copy_demux_fastq,
    demux_source_revision,
    droplet_output_disk_mb,
    init_project_directory,
    print_refresh_result,
    refresh_manifest,
    stage_droplet_reads,
    verify_demux_binary_revision,
    write_demux_completion,
    write_initial_cell_manifest,
    write_run_metadata,
)
from .io import (
    publish,
)
from .maintenance import (
    _ANALYSIS_MASK_INTERVALS,
    append_spikein_contigs,
    begin_reference_build,
    build_downstream_references,
    environment_is_reusable,
    fai_contract_check,
    finalize_release,
    import_reference_binding,
    mask_analysis_reference,
    publish_release,
    record_environment_identity,
    relocate_downstream_bindings,
    seal_downstream_reference_bindings,
    seal_reference_component,
    stage_route_config,
    synthesize_route_config,
    validate_benchmark,
    validate_release,
    verify_route_outputs,
    write_envs_section_spec,
)
from .runtime import (
    read_controller_progress,
)
from .science import (
    _count_bam_pairs,
    call_droplet_cpg,
    emit_rastair_cpg,
    read_pileup_variant_count,
    run_bismark_pair_split,
    run_bsconv_extraction,
    taps_bwa_insert_model,
    write_dupsifter_mqc,
    write_filter_summary,
    write_taps_alignment_qc,
)


# 内部 CLI

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m alopex", description="Alopex 核心模块内部 CLI。")
    commands = parser.add_subparsers(dest="command", required=True)


    init = commands.add_parser("init-project", help="初始化项目骨架（目录/模板/launcher symlink）。")
    init.add_argument("path", nargs="?", default=".", help="项目目录；默认当前目录。")
    init.add_argument("--pipeline-root", required=True, help="Pipeline 安装根，用于模板与 symlink。")
    scalars = commands.add_parser("config-scalars", help="把平台路由所需的 config 标量输出为可 eval 的 shell 赋值。")
    scalars.add_argument("--config", required=True)
    progress = commands.add_parser("controller-progress", help="增量读取 controller 日志进度，不写状态文件。")
    progress.add_argument("--log", required=True)
    progress.add_argument("--offset", type=int, default=0)
    synthesize = commands.add_parser("synthesize-config", help="合成 test 模式 effective config。")
    for field in ('--source', '--destination', '--route'):
        synthesize.add_argument(field, required=True)
    refresh = commands.add_parser("refresh-manifest", help="根据 01_raw/ 刷新 sample manifest。")
    refresh.add_argument("--project", help="项目根目录；默认使用当前目录。")
    refresh.add_argument("--configfile", help="显式指定项目 config；仅接受 <project>/00_config/config.yaml。")
    refresh.add_argument("--dry-run", action="store_true", help="只预览候选 manifest，不写入文件。")


    create = commands.add_parser("snapshot-create", help="创建/复用内容寻址 run snapshot。")
    for field in ('--pipeline-root', '--project', '--config', '--control-env', '--biscuit-env',
                      '--bismark-env', '--rastair-env'):
        create.add_argument(field, required=True)
    create.add_argument("--allow-external-config", action="store_true")
    status = commands.add_parser("project-status", help="输出项目交付状态 JSON。")
    status.add_argument("--project", required=True)
    status.add_argument("--current-snapshot")
    ready = commands.add_parser("delivery-ready", help="写出 delivery_ready 发布事务标记。")
    for field in ('--ready', '--staged-results', '--run-snapshot'):
        ready.add_argument(field, required=True)
    atomic = commands.add_parser("atomic-publish", help="scratch 文件组校验后原子发布。")
    atomic.add_argument("--file", action="append", nargs=2, metavar=("SOURCE", "DESTINATION"),
                        required=True, dest="files")
    atomic.add_argument("--commit-last", default="")
    atomic.add_argument("--bam-destination", default="")
    atomic.add_argument("--samtools", default="samtools")


    pileup_count = commands.add_parser("pileup-variant-count", help="校验当前细胞的 BISCUIT pileup 计数。")
    pileup_count.add_argument("--stats", required=True)
    pileup_count.add_argument("--sample", required=True)

    taps_insert = commands.add_parser("taps-bwa-insert-model", help="推断 TAPS BWA 插入分布并保留长片段。")
    for field in ("r1", "r2", "bwa", "prefix", "scratch"):
        taps_insert.add_argument("--" + field, required=True)
    taps_insert.add_argument("--threads", type=int, required=True)

    taps_qc = commands.add_parser("taps-alignment-qc", help="核算 TAPS BAM 保留与调用集合。")
    for field in ("bam", "samtools", "scratch", "sample", "output"):
        taps_qc.add_argument("--" + field, required=True)
    taps_cpg = commands.add_parser("rastair-cpg", help="转换原生 Rastair 计数并记录对照 QC。")
    for field in ("bed", "reference", "alignment-qc", "output-qc", "bam", "samtools", "scratch", "vcf", "output-snp", "sample"):
        taps_cpg.add_argument("--" + field, required=True)

    write_initial = commands.add_parser("write-cell-manifest", help="写出初始 cell manifest。")
    for field in ('project_manifest', 'barcode_map', 'out_tsv'):
        write_initial.add_argument(field)
    write_initial.add_argument("--species", required=True)
    write_initial.add_argument("--protocol", required=True)
    write_initial.add_argument("--project-sample-id")
    write_initial.add_argument("--called-cells", action="append", type=lambda value: value.split(":", 1), default=[])
    umi = commands.add_parser("call-droplet-cpg", help="按胞嘧啶位点与 UMI 调用分子共识。")
    for field in ("bam", "reference", "coverage", "qc", "mbias", "scratch", "samtools", "sample"):
        umi.add_argument("--" + field, required=True)
    umi.add_argument("--threads", type=int, default=1)
    calling = commands.add_parser("call-droplet-cells", help="按双峰谷底调用 Droplet 细胞。")
    for field in ("counts", "output", "metrics"):
        calling.add_argument("--" + field, required=True)
    calling.add_argument("--demux-binary", required=True)
    calling.add_argument("--threads", type=int, default=1)
    calling.add_argument("--spool", required=True)
    stage_droplet = commands.add_parser("stage-droplet-reads", help="配对解析并提交原始条码暂存与完整计数。")
    for mate in ("r1", "r2"):
        stage_droplet.add_argument("--" + mate, action="append", required=True)
    for field in ("spool", "counts", "structure", "binary", "decoder-python"):
        stage_droplet.add_argument("--" + field, required=True)
    stage_droplet.add_argument("--threads", type=int, required=True)
    write_metadata = commands.add_parser("write-demux-metadata", help="demux 事务内写出 MQC 与 manifest。")
    for arg in ("report_path", "mqc_path", "manifest_path", "dna_dir", "rna_dir",
                "pipeline_mode", "route_id", "downstream_dna", "raw_sample", "cell_manifest_path"):
        write_metadata.add_argument(arg)
    write_metadata.add_argument("--calling-metrics", default="")
    write_metadata.add_argument("--structure-metrics", default="")
    write_completion = commands.add_parser("write-demux-completion", help="demux 事务内写出 completion 提交记录。")
    write_completion.add_argument("out_json")
    for field in ('--run-id', '--project-sample-id', '--manifest', '--report', '--mqc',
                      '--dna-dir', '--rna-dir', '--cell-manifest'):
        write_completion.add_argument(field, required=True)
    collect = commands.add_parser("collect-rna-fastq", help="收集 RNA FASTQ 到保留区。")
    for field in ('--manifest', '--intermediate-dir', '--marker', '--run-id'):
        collect.add_argument(field, required=True)
    write_final = commands.add_parser("write-sample-manifest", help="写出最终 sample manifest。")
    for arg in ("initial_manifest", "intermediate_dir", "staged_results_dir",
                "published_results_dir", "out_tsv", "pair_qc_json", "cutadapt_qc_json"):
        write_final.add_argument(arg)
    write_final.add_argument("--bismark", required=True)
    write_final.add_argument("--workers", required=True, type=int)
    write_final.add_argument("--options-json", required=True)
    write_dupsifter = commands.add_parser("write-dupsifter-qc", help="写出 Dupsifter MultiQC 内容。")
    write_dupsifter.add_argument("stat_path")
    write_dupsifter.add_argument("out_json")
    write_dupsifter.add_argument("--sample", required=True)
    write_dupsifter.add_argument("--metadata", required=True, help="Dupsifter 作业 metadata JSON。")
    extract = commands.add_parser("extract-bsconv", help="bsconv TSV 高 CpH 分类提取。")
    extract.add_argument("--bsconv-tsv", required=True, type=Path)
    extract.add_argument("--threshold", required=True, type=float)
    extract.add_argument("--high-read-names", required=True, type=Path)
    extract.add_argument("--summary-json", required=True, type=Path)
    extract.add_argument("--multiqc-json", type=Path)
    extract.add_argument("--sample-name", required=True)
    extract.add_argument("--excluded-contigs-json", required=True,
                         help="JSON array of canonical reference contig names excluded from assessment.")
    extract.add_argument("--alignment-bam", required=True, type=Path)
    extract.add_argument("--samtools", default="samtools")
    extract.add_argument("--protocol", choices=("cabernet", "srd", "droplet"), required=True)
    extract.add_argument("--high-cph-role", choices=("cdna_contamination", "residual_high_cph"),
                         required=True)
    split = commands.add_parser("split-bismark-pairs", help="Bismark BAM 配对路由与计数。")
    split.add_argument("--excluded-contigs-json", required=True)
    for field in ('--input-bam', '--bio', '--excluded'):
        split.add_argument(field, type=Path, required=True)
    disk = commands.add_parser("droplet-output-disk-mb", help="估算 Droplet 完整文件输出 scratch。")
    disk.add_argument("--spool", required=True, type=Path)
    copy_fastq = commands.add_parser("copy-demux-fastq", help="将完整 FASTQ 从作业 scratch 复制到持久暂存。")
    copy_fastq.add_argument("source", type=Path)
    copy_fastq.add_argument("destination", type=Path)
    copy_fastq.add_argument("--workers", type=int, default=1)
    demux_identity = commands.add_parser("demux-identity", help="计算源码指纹或校验 demux_rs 编译期身份。")
    demux_identity.add_argument("operation", choices=("revision", "verify-binary"))
    demux_identity.add_argument("path")
    demux_identity.add_argument("expected_revision", nargs="?")
    summary = commands.add_parser("write-filter-summary", help="写出 Bismark filter summary（schema v3）。")
    summary.add_argument("--sample", required=True)
    for field in ('--candidate-pairs', '--removed-pairs', '--excluded-pairs'):
        summary.add_argument(field, type=int, required=True)
    summary.add_argument("--processed-bam", type=Path, required=True)
    summary.add_argument("--samtools", required=True)
    summary.add_argument("--percentage", type=float, required=True)
    summary.add_argument("--minimum-count", type=int, required=True)
    summary.add_argument("--protocol", required=True)
    summary.add_argument("--high-cph-role", required=True)
    summary.add_argument("--summary-path", type=Path, required=True)
    summary.add_argument("--multiqc-path", type=Path, required=True)


    bench_validate = commands.add_parser("benchmark-validate", help="校验 Doctor fixture SHA256 与读数守恒。")
    bench_validate.add_argument("--root", required=True)
    bench_validate.add_argument("--manifest", required=True)
    bench_stage = commands.add_parser("benchmark-stage", help="合成沙箱路线 config。")
    for field in ('--project', '--reference', '--route'):
        bench_stage.add_argument(field, required=True)
    bench_verify = commands.add_parser("benchmark-verify", help="校验 sealed 输出守恒。")
    for field in ('--project', '--manifest', '--route'):
        bench_verify.add_argument(field, required=True)
    env_record = commands.add_parser("env-record", help="为一个新 prefix 记录环境身份。")
    env_record.add_argument("prefix")
    env_record.add_argument("--envs", required=True, help="五角色 envs.yaml 路径。")
    env_record.add_argument("--role", required=True)
    env_spec = commands.add_parser("env-spec", help="抽取角色 section 为独立 conda spec 文件。")
    for field in ('--envs', '--role', '--destination'):
        env_spec.add_argument(field, required=True)
    env_finalize = commands.add_parser("env-finalize", help="写出 release manifest。")
    for field in ('release', 'release_id', 'expected_subdir'):
        env_finalize.add_argument(field)
    env_publish = commands.add_parser("env-publish", help="原子切换 conda/current。")
    env_publish.add_argument("current")
    env_publish.add_argument("target")
    env_validate = commands.add_parser("env-validate", help="只读校验 current release。")
    env_validate.add_argument("current")
    env_validate.add_argument("expected_subdir")
    env_validate.add_argument("--envs", required=True)
    env_validate.add_argument("--installed", action="store_true", help="只读核对实际包版本与平台，不导入目标科学模块。")
    env_reusable = commands.add_parser("env-reusable", help="判断既有 prefix 可否复用。")
    env_reusable.add_argument("prefix")
    env_reusable.add_argument("--envs", required=True)
    env_reusable.add_argument("--role", required=True)


    fai_check = commands.add_parser("fai-check", help="FAI 契约校验。")
    fai_check.add_argument("--mode", choices=("source", "assembled"), required=True)
    fai_check.add_argument("--label", required=True)
    fai_check.add_argument("--fai", type=Path, required=True)
    fai_check.add_argument("--fasta", type=Path, required=True)
    suffixes = commands.add_parser("index-suffixes", help="列出 BISCUIT index 后缀清单。")
    suffixes.add_argument("--transient", action="store_true")
    indexes = commands.add_parser("index-paths", help="列出选中 backend/模式的完整 index 文件家族。")
    indexes.add_argument("--reference", required=True)
    indexes.add_argument("--backend", choices=("biscuit", "bismark", "rastair"), required=True)
    indexes.add_argument("--bismark-local-alignment", action="store_true")
    indexes.add_argument("--bismark-root")
    build_downstream = commands.add_parser("build-downstream", help="构建下游 TSS/single-CpG reference。")
    build_downstream.add_argument("--species", required=True, choices=sorted(MANAGED_SPECIES))
    for field in ('--fasta', '--annotation-gtf', '--tss-bed', '--single-cpg-bed'):
        build_downstream.add_argument(field, required=True)
    spikein = commands.add_parser("append-spikein", help="把 lambda/pUC19 追加到 assembled FASTA。")
    for field in ('--lambda-fasta', '--puc19-fasta', '--target-fasta'):
        spikein.add_argument(field, type=Path, required=True)
    mask_par = commands.add_parser("mask-analysis-reference", help="等长屏蔽候选 FASTA 的官方 analysis-set 区间。")
    mask_par.add_argument("--species", choices=sorted(_ANALYSIS_MASK_INTERVALS), required=True)
    mask_par.add_argument("--fasta", type=Path, required=True)
    binding = commands.add_parser("reference-binding", help="Doctor reference 组件的构建输入与内容绑定。")
    binding.add_argument("action", choices=("begin", "seal", "check", "relocate", "import"))
    binding.add_argument("--reference", type=Path, required=True)
    binding.add_argument("--kind", choices=("fasta", "biscuit", "rastair", "bismark-local", "bismark-combined", "downstream"), default="fasta")
    binding.add_argument("--build-input", type=Path)
    binding.add_argument("--root", type=Path)
    binding.add_argument("--annotation-gtf", type=Path)
    binding.add_argument("--tss-bed", type=Path)
    binding.add_argument("--single-cpg-bed", type=Path)
    contigs = commands.add_parser("reference-contigs", help="输出固定assembly的比对contig顺序。")
    contigs.add_argument("--species", required=True, choices=sorted(MANAGED_SPECIES))
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    command = args.command
    try:

        if command == "init-project":
            project = init_project_directory(args.path, args.pipeline_root)
            print(f"Initialized Alopex project: {project}")
            print("Next:")
            print("  1. Put paired FASTQ files under 01_raw/")
            print("  2. ./run_pipeline.sh --refresh-manifest")
            print("  3. ./run_pipeline.sh --dry-run")
            print("  4. ./run_pipeline.sh")
        elif command == "config-scalars":
            config = load_project_config(args.config)
            for var, value in export_config_scalars(config).items():
                print(f"{var}={shlex.quote(value)}")
        elif command == "controller-progress":
            offset, progress = read_controller_progress(args.log, args.offset)
            print(offset, progress)
        elif command == "synthesize-config":
            synthesize_route_config(Path(args.source), Path(args.destination), args.route)
        elif command == "refresh-manifest":
            try:
                project, config_path = _resolve_config(args.project, args.configfile)
                result = refresh_manifest(
                    project_dir=project,
                    config_path=config_path,
                    dry_run=bool(args.dry_run),
                )
            except Exception as exc:
                print(f"ERROR: {exc}", file=sys.stderr)
                return 2
            return print_refresh_result(args, result)

        elif command == "snapshot-create":
            print(create_run_snapshot(
                pipeline_root=args.pipeline_root,
                project_dir=args.project,
                config_path=args.config,
                control_env=args.control_env,
                biscuit_env=args.biscuit_env,
                bismark_env=args.bismark_env,
                rastair_env=args.rastair_env,
                allow_external_config=args.allow_external_config,
            ))
        elif command == "project-status":
            result = project_status(
                project_dir=args.project,
                current_snapshot=args.current_snapshot,
            )
            print(json.dumps(result, ensure_ascii=False, sort_keys=True))
        elif command == "delivery-ready":
            run_snapshot = str(args.run_snapshot).strip()
            if not run_snapshot:
                raise ValueError(
                    "publishing requires a launcher-managed run snapshot. "
                    "Direct Snakemake runs stop at 02_work/results_stage; rerun via "
                    "core/run_pipeline.sh to publish."
                )
            write_delivery_ready(
                args.ready,
                staged_results_dir=args.staged_results,
                run_snapshot_path=run_snapshot,
            )
        elif command == "atomic-publish":
            pairs = [(Path(source), Path(destination)) for source, destination in args.files]
            publish(
                pairs,
                commit_last=Path(args.commit_last) if args.commit_last else None,
                bam_destination=Path(args.bam_destination) if args.bam_destination else None,
                samtools=args.samtools,
            )

        elif command == "pileup-variant-count":
            print(read_pileup_variant_count(args.stats, args.sample))
        elif command == "call-droplet-cpg":
            call_droplet_cpg(args.bam, args.reference, args.coverage, args.qc, args.mbias, sample=args.sample, scratch=args.scratch, samtools=args.samtools, threads=args.threads)
        elif command == "call-droplet-cells":
            call_droplet_cells(args.counts, args.output, args.metrics, demux_binary=args.demux_binary, threads=args.threads, spool=args.spool)
        elif command == "stage-droplet-reads":
            stage_droplet_reads(args.r1, args.r2, args.spool, args.counts, args.structure,
                                binary=args.binary, decoder_python=args.decoder_python, threads=args.threads)
        elif command == "write-cell-manifest":
            write_initial_cell_manifest(
                args.project_manifest, args.barcode_map, args.out_tsv,
                species=args.species, protocol=args.protocol, called_cells=args.called_cells,
                project_sample_id=args.project_sample_id,
            )
        elif command == "write-demux-metadata":
            write_run_metadata(
                args.report_path, args.mqc_path, args.manifest_path, args.dna_dir,
                args.rna_dir, args.pipeline_mode, args.route_id, args.downstream_dna,
                args.raw_sample, args.cell_manifest_path,
                calling_metrics=args.calling_metrics, structure_metrics=args.structure_metrics,
            )
        elif command == "write-demux-completion":
            write_demux_completion(
                args.out_json,
                run_id=args.run_id,
                project_sample_id=args.project_sample_id,
                manifest_path=args.manifest,
                report_path=args.report,
                mqc_path=args.mqc,
                dna_dir=args.dna_dir,
                rna_dir=args.rna_dir,
                cell_manifest_path=args.cell_manifest,
            )
        elif command == "collect-rna-fastq":
            collect_rna_fastq(args.manifest, args.intermediate_dir, args.marker, args.run_id)
        elif command == "write-sample-manifest":
            write_final_sample_manifest(
                args.initial_manifest, args.intermediate_dir, args.staged_results_dir,
                args.out_tsv,
                options=json.loads(args.options_json),
                published_results_dir=args.published_results_dir,
                bismark=args.bismark,
                pair_qc_multiqc_json=args.pair_qc_json,
                cutadapt_qc_multiqc_json=args.cutadapt_qc_json,
                workers=args.workers,
            )
        elif command == "taps-bwa-insert-model":
            print(taps_bwa_insert_model(args.r1, args.r2, args.bwa, args.prefix, args.scratch, args.threads))
        elif command == "taps-alignment-qc":
            write_taps_alignment_qc(args.bam, args.samtools, args.scratch, args.sample, args.output)
        elif command == "rastair-cpg":
            emit_rastair_cpg(args.bed, args.reference, args.alignment_qc, args.output_qc, sys.stdout,
                             bam=args.bam, samtools=args.samtools, scratch=args.scratch,
                             vcf=args.vcf, output_snp=args.output_snp, sample=args.sample)
        elif command == "write-dupsifter-qc":
            write_dupsifter_mqc(
                args.stat_path, args.out_json, args.sample,
                metadata=json.loads(args.metadata),
            )
        elif command == "extract-bsconv":
            run_bsconv_extraction(args)
        elif command == "split-bismark-pairs":
            return run_bismark_pair_split(args)
        elif command == "droplet-output-disk-mb":
            print(droplet_output_disk_mb(args.spool))
        elif command == "copy-demux-fastq":
            print(json.dumps(copy_demux_fastq(args.source, args.destination, args.workers)))
        elif command == "demux-identity":
            if args.operation == "revision":
                if args.expected_revision is not None:
                    raise ValueError("demux-identity revision accepts only the source path")
                print(demux_source_revision(args.path))
            else:
                if args.expected_revision is None:
                    raise ValueError("demux-identity verify-binary requires the expected revision")
                verify_demux_binary_revision(args.path, args.expected_revision)
        elif command == "write-filter-summary":
            write_filter_summary(
                sample=args.sample,
                candidate_pairs=args.candidate_pairs,
                removed_pairs=args.removed_pairs,
                excluded_pairs=args.excluded_pairs,
                retained_pairs=_count_bam_pairs(args.processed_bam, args.samtools),
                percentage=args.percentage,
                minimum_count=args.minimum_count,
                protocol=args.protocol,
                high_cph_role=args.high_cph_role,
                summary=args.summary_path,
                multiqc=args.multiqc_path,
            )

        elif command == "benchmark-validate":
            validate_benchmark(Path(args.root), Path(args.manifest))
        elif command == "benchmark-stage":
            stage_route_config(Path(args.project), Path(args.reference), args.route)
        elif command == "benchmark-verify":
            verify_route_outputs(Path(args.project), Path(args.manifest), args.route)
        elif command == "env-record":
            record_environment_identity(
                prefix=args.prefix,
                role=args.role,
                spec_name=envs_spec_logical_name(args.role),
                spec_sha256=envs_section_sha256(args.envs, args.role),
            )
        elif command == "env-spec":
            write_envs_section_spec(args.envs, args.role, args.destination)
        elif command == "env-finalize":
            finalize_release(args.release, args.release_id, args.expected_subdir)
        elif command == "env-publish":
            publish_release(args.current, args.target)
        elif command == "env-validate":
            validate_release(args.current, args.expected_subdir, args.envs, installed=args.installed)
        elif command == "env-reusable":
            environment_is_reusable(args.prefix, args.role, args.envs)

        elif command == "reference-binding":
            if args.action == "begin":
                if args.build_input is None:
                    raise ValueError("--build-input is required")
                begin_reference_build(args.reference, args.build_input, args.annotation_gtf)
            elif args.kind == "downstream":
                if args.tss_bed is None or args.single_cpg_bed is None:
                    raise ValueError("downstream bindings require TSS and single-CpG paths")
                if args.action == "relocate":
                    relocate_downstream_bindings(args.reference, args.tss_bed, args.single_cpg_bed)
                elif args.action == "import":
                    for path, kind in ((args.tss_bed, "tss"), (args.single_cpg_bed, "single-cpg")):
                        import_reference_binding(args.reference, kind, artifact=path, annotation=args.annotation_gtf)
                elif args.action == "seal":
                    seal_downstream_reference_bindings(args.reference, args.annotation_gtf,
                        args.tss_bed, args.single_cpg_bed, args.build_input)
                else:
                    validate_downstream_reference(args.reference, args.tss_bed, "tss", args.annotation_gtf)
                    validate_downstream_reference(args.reference, args.single_cpg_bed, "single-cpg", args.annotation_gtf)
            elif args.action == "import":
                import_reference_binding(args.reference, args.kind, root=args.root)
            elif args.action == "seal":
                seal_reference_component(args.reference, args.kind, build_input=args.build_input, root=args.root)
            else:
                validate_reference_component(args.reference, args.kind, root=args.root)
        elif command == "reference-contigs":
            print("\n".join(name for name, _length in managed_alignment_contigs(args.species) if name not in SPIKE_CONTIGS))
        elif command == "fai-check":
            return 0 if fai_contract_check(args.mode, args.label, args.fai, args.fasta) else 1
        elif command == "index-suffixes":
            suffixes_list = BISCUIT_TRANSIENT_SUFFIXES if args.transient else BISCUIT_INDEX_SUFFIXES
            print("\n".join(suffixes_list))
        elif command == "index-paths":
            print("\n".join(map(str, backend_active_index_paths(
                args.reference, args.backend,
                bismark_local_alignment=args.bismark_local_alignment,
                bismark_root=args.bismark_root,
            ))))
        elif command == "append-spikein":
            append_spikein_contigs(args.lambda_fasta, args.puc19_fasta, args.target_fasta)
        elif command == "mask-analysis-reference":
            mask_analysis_reference(args.species, args.fasta)
        elif command == "build-downstream":
            payload = build_downstream_references(
                species=args.species,
                fasta=args.fasta,
                annotation_gtf=args.annotation_gtf,
                tss_bed=args.tss_bed,
                single_cpg_bed=args.single_cpg_bed,
            )
            print(json.dumps(payload, ensure_ascii=False, sort_keys=True))
        else:
            raise ValueError(f"unsupported alopex command: {command}")
    except (
        ValueError,
        RuntimeError,
        OSError,
        subprocess.SubprocessError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

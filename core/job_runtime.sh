#!/usr/bin/env bash


if [[ -n "${DNA_PIPELINE_ROOT:-}" ]]; then
    export PYTHONPATH="$DNA_PIPELINE_ROOT/core${PYTHONPATH:+:$PYTHONPATH}"
fi

# 首次尝试即建立日志硬链接，使 Snakemake 重试清理主日志后仍保留失败证据。
dna_pipeline_init_attempt_logs() {
    local attempt="$1"
    shift
    local live_log runtime_host

    [[ "$attempt" =~ ^[12]$ ]] || {
        echo "ERROR: unsupported attempt number: $attempt" >&2
        return 2
    }
    (( $# > 0 )) || {
        echo "ERROR: at least one live log path is required" >&2
        return 2
    }
    runtime_host="$(hostname)"
    for live_log in "$@"; do
        mkdir -p "$(dirname "$live_log")"
        rm -f -- "$live_log"
        if [[ "$attempt" == "1" ]]; then
            rm -f -- "${live_log}.attempt1"
            : > "${live_log}.attempt1"
            ln -- "${live_log}.attempt1" "$live_log"
        else
            : > "$live_log"
        fi
        printf '=== resource attempt %s ===\n' "$attempt" >> "$live_log"
        printf 'runtime_host=%s slurm_job_id=%s slurm_partition=%s slurm_qos=%s slurm_cpus_per_task=%s slurm_mem_per_node=%s\n' \
            "$runtime_host" "${SLURM_JOB_ID:-NA}" "${SLURM_JOB_PARTITION:-NA}" \
            "${SLURM_JOB_QOS:-NA}" "${SLURM_CPUS_PER_TASK:-NA}" \
            "${SLURM_MEM_PER_NODE:-NA}" >> "$live_log"
    done
}

# 在有空间余量的首选候选盘上创建本 job 专属 scratch 目录并导出为 TMPDIR；请求 110% 余量且路径锁定在 namespace 内。
dna_pipeline_allocate_scratch() {
    local configured_tmp="$1"
    local namespace="$2"
    local group="$3"
    local cell="$4"
    local disk_mb="$5"
    local project_root="${6:-${DNA_PIPELINE_PROJECT_DIR:-}}"
    local node_override="${DNA_PIPELINE_NODE_TMPDIR:-}"
    local project_tmp="${project_root:+$project_root/05_tmp}"
    local base_tmp="" candidate available_kb required_kb scratch_parent scratch_suffix token

    for token in "$namespace" "$group" "$cell"; do
        [[ -n "$token" && "$token" != . && "$token" != .. && "$token" != *[!A-Za-z0-9_.-]* ]] || {
            echo "ERROR: unsafe scratch path component: $token" >&2
            return 2
        }
    done
    [[ "$disk_mb" =~ ^[0-9]+$ ]] || {
        echo "ERROR: invalid scratch request: $disk_mb MiB" >&2
        return 2
    }
    [[ -n "$project_root" && "$project_root" == /* && "$project_root" != "/" ]] || {
        echo "ERROR: an absolute non-root project path is required for scratch fallback" >&2
        return 2
    }
    mkdir -p "$project_tmp"
    required_kb=$(((disk_mb * 1024 * 110 + 99) / 100))

    for candidate in "${SLURM_TMPDIR:-}" "$node_override" "$configured_tmp" "$project_tmp"; do
        [[ -n "$candidate" && "$candidate" == /* && "$candidate" != "/" ]] || continue
        [[ -d "$candidate" && -w "$candidate" ]] || continue
        candidate="$(cd "$candidate" && pwd -P)"
        available_kb="$(df -Pk "$candidate" | awk 'END{print $4}')"
        [[ "$available_kb" =~ ^[0-9]+$ ]] || continue
        if (( available_kb < required_kb )); then
            echo "WARNING: scratch candidate has insufficient free space: path=$candidate available_kb=$available_kb required_with_margin_kb=$required_kb" >&2
            continue
        fi
        base_tmp="$candidate"
        scratch_suffix="${SLURM_JOB_ID:-local}.$$"
        [[ "$scratch_suffix" != *[!A-Za-z0-9_.-]* ]] || scratch_suffix="local.$$"
        scratch_parent="$candidate"
        for token in dna_pipeline "$namespace" "$group"; do
            scratch_parent="$scratch_parent/$token"
            [[ ! -L "$scratch_parent" && ( ! -e "$scratch_parent" || -d "$scratch_parent" ) ]] || {
                echo "ERROR: scratch parent must be a real directory: $scratch_parent" >&2
                return 2
            }
            mkdir -p "$scratch_parent" || return $?
        done
        DNA_PIPELINE_SCRATCH_DIR="$(mktemp -d "$scratch_parent/$cell.$scratch_suffix.XXXXXX")" || return $?
        break
    done
    [[ -n "$base_tmp" ]] || {
        echo "ERROR: no writable scratch candidate has the requested space plus 10% margin" >&2
        return 1
    }

    case "$DNA_PIPELINE_SCRATCH_DIR" in
        "$base_tmp"/dna_pipeline/"$namespace"/"$group"/"$cell".*) ;;
        *) echo "ERROR: resolved scratch path escaped its namespace" >&2; return 2 ;;
    esac
    DNA_PIPELINE_SCRATCH_BASE="$base_tmp"
    TMPDIR="$DNA_PIPELINE_SCRATCH_DIR"
    export DNA_PIPELINE_SCRATCH_BASE DNA_PIPELINE_SCRATCH_DIR TMPDIR
}

# 安装 EXIT trap：job 以任何方式退出（含失败）时自动清理 scratch。
dna_pipeline_install_scratch_cleanup() {
    trap dna_pipeline_cleanup_scratch_on_exit EXIT
}

# 删除当前 scratch 目录；路径逃逸本 namespace 时拒绝执行并报错。
dna_pipeline_cleanup_scratch() {
    local base="${DNA_PIPELINE_SCRATCH_BASE:-}"
    local scratch="${DNA_PIPELINE_SCRATCH_DIR:-}"
    [[ -n "$scratch" ]] || return 0
    case "$scratch" in
        "$base"/dna_pipeline/*/*/*) rm -rf -- "$scratch" ;;
        *)
            echo "ERROR: refusing to clean unsafe scratch path: $scratch" >&2
            return 2
            ;;
    esac
}

# trap 清理入口：保留原始退出码，执行 scratch 清理后按原码退出。
dna_pipeline_cleanup_scratch_on_exit() {
    local exit_status=$?
    trap - EXIT
    dna_pipeline_cleanup_scratch || true
    exit "$exit_status"
}

# 成功路径的显式释放：取消 EXIT trap 并立即清理 scratch。
dna_pipeline_release_scratch() {
    trap - EXIT
    dna_pipeline_cleanup_scratch
}

# 向调用方重定向的日志输出 hostname、SLURM job 与 scratch 路径/挂载及容量信息。
dna_pipeline_log_scratch() {
    printf 'INFO: hostname=%s slurm_job_id=%s\n' "$(hostname)" "${SLURM_JOB_ID:-NA}"
    printf 'INFO: scratch_path=%s scratch_mount=%s\n' \
        "$DNA_PIPELINE_SCRATCH_DIR" \
        "$(df -P "$DNA_PIPELINE_SCRATCH_BASE" | awk 'END{print $6}')"
    df -h "$DNA_PIPELINE_SCRATCH_BASE"
}

# --recursive 递归搜索；--exclude-glob 可多次给出，命中文件名即从候选中剔除。
dna_pipeline_find_single_output() {
    local directory="$1"
    shift
    local recursive=0
    local -a exclude_globs=()
    while [[ "${1:-}" == --* ]]; do
        case "$1" in
            --recursive) recursive=1 ;;
            --exclude-glob)
                [[ -n "${2:-}" ]] || { echo "ERROR: --exclude-glob requires a pattern" >&2; return 2; }
                exclude_globs+=("$2")
                shift
                ;;
            *) echo "ERROR: unsupported option: $1" >&2; return 2 ;;
        esac
        shift
    done
    local -a patterns=("$@")
    local -a matches=()
    local pattern path glob excluded
    shopt -s nullglob
    for pattern in "${patterns[@]}"; do
        if (( recursive == 1 )); then
            while IFS= read -r -d '' path; do matches+=("$path"); done < <(find "$directory" -type f -name "$pattern" -print0)
        else
            while IFS= read -r -d '' path; do matches+=("$path"); done < <(find "$directory" -maxdepth 1 -type f -name "$pattern" -print0)
        fi
    done
    shopt -u nullglob
    if (( ${#exclude_globs[@]} > 0 )); then
        local -a kept=()
        for path in ${matches[@]+"${matches[@]}"}; do
            excluded=0
            for glob in "${exclude_globs[@]}"; do
                if [[ "${path##*/}" == $glob ]]; then
                    excluded=1
                    break
                fi
            done
            (( excluded == 0 )) && kept+=("$path")
        done
        matches=(${kept[@]+"${kept[@]}"})
    fi
    local -a unique=()
    for path in ${matches[@]+"${matches[@]}"}; do
        local seen=0 item
        for item in ${unique[@]+"${unique[@]}"}; do [[ "$item" == "$path" ]] && seen=1 && break; done
        (( seen == 0 )) && unique+=("$path")
    done
    local unique_count=${#unique[@]}
    if (( unique_count != 1 )); then
        printf 'ERROR: expected exactly one matching file in %s; observed=%s\n' \
            "$directory" "${unique[*]:-none}" >&2
        return 2
    fi
    printf '%s\n' "${unique[0]}"
}


# 把 mapped/variant 计数序列化为 final manifest 直接消费的稳定 pileup stats JSON。
dna_pipeline_write_pileup_stats() {
    local output="$1" sample="$2" mapped="$3" variants="$4" empty=false
    [[ "$mapped" =~ ^[0-9]+$ && "$variants" =~ ^[0-9]+$ ]] || {
        echo "ERROR: pileup stats counts must be non-negative integers" >&2
        return 2
    }
    (( variants == 0 )) && empty=true
    mkdir -p "$(dirname "$output")"
    printf '{"schema_version":1,"sample":"%s","mapped_reads":%s,"variant_records":%s,"empty":%s}\n' \
        "$sample" "$mapped" "$variants" "$empty" > "$output"
}

# 读取 pileup stats 的指定计数字段；缺失或无效字段直接失败。
dna_pipeline_pileup_stat() {
    local path="$1" field="$2" value
    case "$field" in
        mapped_reads|variant_records)
            value="$(sed -n "s/.*\\\"${field}\\\":[[:space:]]*\\([0-9][0-9]*\\).*/\\1/p" "$path")"
            ;;
        empty)
            value="$(sed -n 's/.*"empty":[[:space:]]*\(true\|false\).*/\1/p' "$path")"
            ;;
        *) echo "ERROR: unsupported pileup stats field: $field" >&2; return 2 ;;
    esac
    [[ -n "$value" ]] || { echo "ERROR: malformed pileup stats: $path" >&2; return 1; }
    printf '%s\n' "$value"
}

# JSON 由 cutadapt 自身写出并以退出码保证；发布侧 atomic-publish 再做 JSON 解析校验。
dna_pipeline_run_cutadapt() {
    local cutadapt="$1" gzip_bin="$2" threads="$3" quality="$4" min_length="$5"
    local adapter_r1="$6" adapter_r2="$7" r1_in="$8" r2_in="$9"
    local r1_out="${10}" r2_out="${11}" json_out="${12}" log_file="${13}"
    "$cutadapt" \
        -a "$adapter_r1" -A "$adapter_r2" \
        --error-rate 0.2 \
        -q "$quality" -m "$min_length" -j "$threads" \
        -o "$r1_out" -p "$r2_out" --json "$json_out" \
        "$r1_in" "$r2_in" >> "$log_file" 2>&1
    "$gzip_bin" -t "$r1_out" "$r2_out" >> "$log_file" 2>&1
}

# 运行 high-CpH 分类：零读 BAM 直接产出空集合，否则 biscuit bsconv -p 流式进入解析器；两级 PIPESTATUS 任一非零即失败。
dna_pipeline_run_bsconv_classification() {
    local samtools="$1" biscuit="$2" ref="$3" bam="$4" log_file="$5"
    shift 5
    local pipeline_src="${DNA_PIPELINE_ROOT:?}/core"
    local -a parser_args=("$@")
    local total_reads
    total_reads="$("$samtools" idxstats "$bam" | awk '{n += $3 + $4} END {print n + 0}')"
    if [[ "$total_reads" -eq 0 ]]; then
        echo "WARNING: BAM contains zero reads; high-CpH classification is empty." >> "$log_file"
        PYTHONPATH="$pipeline_src" python -m dna_pipeline extract-bsconv "${parser_args[@]}" < /dev/null >> "$log_file" 2>&1
    else
        set +e
        "$biscuit" bsconv -p "$ref" "$bam" 2>> "$log_file" | \
            PYTHONPATH="$pipeline_src" python -m dna_pipeline extract-bsconv "${parser_args[@]}" >> "$log_file" 2>&1
        local -a status=("${PIPESTATUS[@]}")
        set -e
        if (( ${status[0]:-99} != 0 || ${status[1]:-99} != 0 )); then
            echo "ERROR: high-CpH stream failed: biscuit=${status[0]:-99} parser=${status[1]:-99}" >> "$log_file"
            exit 1
        fi
    fi
}

# 把后端计数流统一写成四列 canonical CpG；压缩固定单线程，避免嵌套并行。
dna_pipeline_emit_canonical_cpg() {
    local zstd="${1:?}" pos_expr="${2:?}" out_tmp="${3:?}"
    {
        printf 'chrom\tpos\tmethyl\tunmethyl\n'
        awk -v OFF="$pos_expr" 'BEGIN{OFS="\t"} NF>=6 {
            pos=$2+OFF; m=$5+0; u=$6+0;
            if(pos>=0 && m+u>0) print $1,pos,m,u;
        }'
    } | "$zstd" -q -T1 -c > "$out_tmp"
}

# 对 BAM 坐标排序，并通过 atomic-publish 以 BAM commit-last 原子提交。
dna_pipeline_coordinate_sort_and_stage() {
    local samtools="${1:?}" extra="${2:?}" src="${3:?}"
    local tmp_bam="${4:?}" tmp_bai="${5:?}" final_bam="${6:?}" final_bai="${7:?}"
    mkdir -p "$(dirname "$final_bam")"
    "$samtools" sort -@ "$extra" -o "$tmp_bam" "$src"
    "$samtools" index -@ "$extra" "$tmp_bam"
    "$samtools" quickcheck -v "$tmp_bam"
    python -m dna_pipeline atomic-publish \
        --samtools "$samtools" \
        --bam-destination "$final_bam" \
        --commit-last "$final_bam" \
        --file "$tmp_bai" "$final_bai" \
        --file "$tmp_bam" "$final_bam"
}

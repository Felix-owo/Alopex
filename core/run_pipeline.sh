#!/usr/bin/env bash
set -euo pipefail

# 按脚本真实路径定位 Pipeline，支持 launcher 符号链接。
dna_pipeline_resolve_launcher_source() {
    local source="${BASH_SOURCE[0]}" directory target
    while [[ -L "$source" ]]; do
        directory="$(cd -P "$(dirname "$source")" && pwd)"
        target="$(readlink "$source")"
        if [[ "$target" == /* ]]; then
            source="$target"
        else
            source="$directory/$target"
        fi
    done
    directory="$(cd -P "$(dirname "$source")" && pwd)"
    printf '%s/%s\n' "$directory" "$(basename "$source")"
}
DNA_PIPELINE_LAUNCHER_SOURCE="$(dna_pipeline_resolve_launcher_source)"
DNA_PIPELINE_ROOT="$(cd "$(dirname "$DNA_PIPELINE_LAUNCHER_SOURCE")/.." && pwd -P)"
export DNA_PIPELINE_ROOT

source "$DNA_PIPELINE_ROOT/core/doctor.sh"

unset DNA_PIPELINE_RUN_SNAPSHOT DNA_PIPELINE_DELIVERY_ID
DNA_PIPELINE_DATA_MODE="${DNA_PIPELINE_DATA_MODE:-production}"
DNA_PIPELINE_TEST_ROUTE="${DNA_PIPELINE_TEST_ROUTE:-auto}"
DNA_PIPELINE_DOCTOR_LOCAL_CORES="${DNA_PIPELINE_DOCTOR_LOCAL_CORES:-}"
DNA_PIPELINE_INIT_PROJECT=""
DNA_PIPELINE_REFRESH_MANIFEST=0
DNA_PIPELINE_READ_ONLY=0

usage() {
    cat <<'EOF'
Alopex 项目启动器

常规用户流程（在项目目录内运行）：
  ./run_pipeline.sh --init-project
  ./run_pipeline.sh --refresh-manifest
  ./run_pipeline.sh --dry-run
  ./run_pipeline.sh

项目发现：
  当前目录包含 00_config/config.yaml 时直接使用；否则自动向上层目录查找。

选项：
  --init-project [PROJECT]  初始化项目骨架（独立维护动作，完成后立即退出）。
  --refresh-manifest        重扫 01_raw/ 并刷新 sample manifest（独立维护动作）。
  --dry-run, -n             只读检查：不写任何项目状态。
  -h, --help                显示本帮助。

平台路由：
  macOS                         -> Snakemake local 执行器
  Linux + Slurm login 节点      -> 当前 shell 作为 controller，Snakemake Slurm executor 提交 rule jobs
  Linux + 已有 SLURM allocation -> 当前计算节点直接作为 controller，继续通过 Slurm executor 提交 rule jobs

启动器不接受其他选项或 Snakemake target；未知参数一律拒绝。
EOF
}

while (( $# > 0 )); do
    case "$1" in
        --help|-h)
            usage; exit 0 ;;
        --init-project)
            if (( $# >= 2 )) && [[ "$2" != -* ]]; then
                [[ -n "$2" ]] || { echo "ERROR: --init-project requires a non-empty path." >&2; exit 2; }
                DNA_PIPELINE_INIT_PROJECT="$2"; shift 2
            else
                DNA_PIPELINE_INIT_PROJECT="$(pwd -P)"; shift
            fi ;;
        --init-project=*)
            DNA_PIPELINE_INIT_PROJECT="${1#--init-project=}"
            [[ -n "$DNA_PIPELINE_INIT_PROJECT" ]] || { echo "ERROR: --init-project requires a non-empty path." >&2; exit 2; }
            shift ;;
        --refresh-manifest)
            DNA_PIPELINE_REFRESH_MANIFEST=1; shift ;;
        --dry-run|-n)
            DNA_PIPELINE_READ_ONLY=1; shift ;;
        -*)
            echo "ERROR: unsupported launcher argument: $1 (see --help)" >&2; exit 2 ;;
        *)
            echo "ERROR: launcher accepts no Snakemake targets or extra arguments: $1 (see --help)" >&2; exit 2 ;;
    esac
done

case "$DNA_PIPELINE_DATA_MODE" in production|test) ;; *) echo "ERROR: DNA_PIPELINE_DATA_MODE must be production or test." >&2; exit 2;; esac
case "$DNA_PIPELINE_TEST_ROUTE" in auto|hg38-cabernet-bismark|mm10-srd-biscuit|hg38-taps-rastair) ;; *) echo "ERROR: unsupported DNA_PIPELINE_TEST_ROUTE: $DNA_PIPELINE_TEST_ROUTE" >&2; exit 2;; esac
if [[ "$DNA_PIPELINE_DATA_MODE" == production && "$DNA_PIPELINE_TEST_ROUTE" != auto ]]; then
    echo "ERROR: DNA_PIPELINE_TEST_ROUTE is only valid when DNA_PIPELINE_DATA_MODE=test." >&2; exit 2
fi
if [[ -n "$DNA_PIPELINE_DOCTOR_LOCAL_CORES" ]]; then
    [[ "$DNA_PIPELINE_DATA_MODE" == test && "$DNA_PIPELINE_TEST_ROUTE" != auto ]] || {
        echo "ERROR: DNA_PIPELINE_DOCTOR_LOCAL_CORES is only valid for an internal Doctor test route." >&2; exit 2
    }
    [[ "$DNA_PIPELINE_DOCTOR_LOCAL_CORES" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: DNA_PIPELINE_DOCTOR_LOCAL_CORES must be positive." >&2; exit 2
    }
fi


if (( DNA_PIPELINE_REFRESH_MANIFEST == 1 )); then
    (( DNA_PIPELINE_READ_ONLY == 0 )) || {
        echo "ERROR: --refresh-manifest is a standalone maintenance action and cannot be combined with --dry-run." >&2
        echo "Run './run_pipeline.sh --refresh-manifest' first, then run './run_pipeline.sh --dry-run' separately." >&2
        exit 2
    }
fi

DNA_PIPELINE_LAUNCH_DIR="$(pwd -P)"

DNA_PIPELINE_CONTROL_ENV="$DNA_PIPELINE_ROOT/conda/current/control"

if [[ -n "$DNA_PIPELINE_INIT_PROJECT" ]]; then
    (( DNA_PIPELINE_REFRESH_MANIFEST == 0 )) || { echo "ERROR: --init-project cannot be combined with --refresh-manifest." >&2; exit 2; }
    (( DNA_PIPELINE_READ_ONLY == 0 )) || { echo "ERROR: --init-project cannot be combined with --dry-run." >&2; exit 2; }
    if [[ -x "$DNA_PIPELINE_CONTROL_ENV/bin/python" ]]; then
        DNA_PIPELINE_INIT_PYTHON="$DNA_PIPELINE_CONTROL_ENV/bin/python"
    else
        DNA_PIPELINE_INIT_PYTHON="$(command -v python3 || command -v python)"
    fi
    PYTHONPATH="$DNA_PIPELINE_ROOT/core" "$DNA_PIPELINE_INIT_PYTHON" -m dna_pipeline init-project \
        "$DNA_PIPELINE_INIT_PROJECT" --pipeline-root "$DNA_PIPELINE_ROOT"
    exit $?
fi
DNA_PIPELINE_BISCUIT_ENV="$DNA_PIPELINE_ROOT/conda/current/biscuit"
DNA_PIPELINE_BISMARK_ENV="$DNA_PIPELINE_ROOT/conda/current/bismark"
DNA_PIPELINE_RASTAIR_ENV="$DNA_PIPELINE_ROOT/conda/current/rastair"
DNA_PIPELINE_RULE_ENV=""

while IFS= read -r variable_name; do
    case "$variable_name" in
        SBATCH_*|SNAKEMAKE_*)
            unset "$variable_name"
            ;;
    esac
done < <(compgen -e)
unset DNA_PIPELINE_NODE_TMPDIR

require_environment_tools() {
    local prefix="$1" role="$2" relative
    shift 2
    for relative in "$@"; do
        [[ -x "$prefix/$relative" ]] || {
            echo "ERROR: $role environment is missing executable: $prefix/$relative" >&2
            echo "Repair it with: $DNA_PIPELINE_ROOT/core/doctor.sh" >&2
            exit 127
        }
    done
}

require_environment_tools "$DNA_PIPELINE_CONTROL_ENV" control "${CONTROL_RUNTIME_REQUIRED[@]}"
DNA_PIPELINE_CONTROL_ENV="$(cd "$DNA_PIPELINE_CONTROL_ENV" && pwd -P)"
DNA_PIPELINE_PYTHON="$DNA_PIPELINE_CONTROL_ENV/bin/python"
DNA_PIPELINE_SNAKEMAKE="$DNA_PIPELINE_CONTROL_ENV/bin/snakemake"
DNA_PIPELINE_CONDA="$DNA_PIPELINE_CONTROL_ENV/bin/conda"
export DNA_PIPELINE_CONTROL_ENV

# 从当前目录逐级向上查找最近的 <project>/00_config/config.yaml。
discover_project_config() {
    local cursor="$1" candidate parent
    cursor="$(cd "$cursor" && pwd -P)"
    while :; do
        candidate="$cursor/00_config/config.yaml"
        if [[ -f "$candidate" ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
        parent="$(dirname "$cursor")"
        [[ "$parent" != "$cursor" ]] || break
        cursor="$parent"
    done
    return 1
}

DNA_PIPELINE_CONFIG="$(discover_project_config "$DNA_PIPELINE_LAUNCH_DIR" || true)"
[[ -f "$DNA_PIPELINE_CONFIG" ]] || {
    echo "ERROR: no Alopex project found in the current directory or its parents." >&2
    echo "Initialize the current directory with: /path/to/Alopex_dev/core/run_pipeline.sh --init-project" >&2
    exit 2
}
DNA_PIPELINE_PROJECT_DIR="$(cd "$(dirname "$DNA_PIPELINE_CONFIG")/.." && pwd -P)"
export DNA_PIPELINE_CONFIG DNA_PIPELINE_PROJECT_DIR

# 清理 launcher 临时合成目录。
cleanup_launcher_tmp() {
    if [[ -n "${DNA_PIPELINE_RUNTIME_TMP:-}" && -d "${DNA_PIPELINE_RUNTIME_TMP:-}" ]]; then
        case "$(basename "$DNA_PIPELINE_RUNTIME_TMP")" in
            dna-pipeline-launch.*) rm -rf -- "$DNA_PIPELINE_RUNTIME_TMP" ;;
            *) printf 'WARNING: refusing to remove unexpected launcher temp path: %s\n' "$DNA_PIPELINE_RUNTIME_TMP" >&2 ;;
        esac
    fi
    DNA_PIPELINE_RUNTIME_TMP=""
}

if (( DNA_PIPELINE_REFRESH_MANIFEST == 1 )); then
    unset PYTHONHOME
    export PYTHONNOUSERSITE=1
    export PYTHONPATH="$DNA_PIPELINE_ROOT/core"
    "$DNA_PIPELINE_PYTHON" -m dna_pipeline refresh-manifest --project "$DNA_PIPELINE_PROJECT_DIR"
    exit $?
fi

project_hash8() {
    "$DNA_PIPELINE_PYTHON" -c 'import hashlib,sys;print(hashlib.sha256(sys.argv[1].encode()).hexdigest()[:8])' "$DNA_PIPELINE_PROJECT_DIR"
}

pipeline_interrupt_tracked_run() {
    local exit_code="$1"
    if [[ "$DNA_PIPELINE_ACTIVE_SNAKEMAKE_PID" =~ ^[0-9]+$ ]]; then
        txn_terminate_pid_bounded \
            "$DNA_PIPELINE_ACTIVE_SNAKEMAKE_PID" \
            "Snakemake controller" \
            "${DNA_PIPELINE_TERMINATION_GRACE_SECONDS:-15}"
        DNA_PIPELINE_ACTIVE_SNAKEMAKE_PID=""
    fi
    exit "$exit_code"
}
# 受控执行 Snakemake：接管中断信号；完整输出进 controller log，进度以 Snakemake 原生输出为准。
run_snakemake_tracked() {
    local run_id="$1" rc=0 controller_log
    local progress_interval last done_n total_n pct_n
    shift
    local arg is_dry_run=0
    for arg in "$@"; do [[ "$arg" == "--dry-run" ]] && is_dry_run=1; done
    if (( DNA_PIPELINE_READ_ONLY == 1 || is_dry_run == 1 )); then
        "$@"
        return $?
    fi
    trap 'pipeline_interrupt_tracked_run 130' INT
    trap 'pipeline_interrupt_tracked_run 143' TERM
    trap 'pipeline_interrupt_tracked_run 129' HUP
    set +e
    controller_log="$DNA_PIPELINE_PROJECT_DIR/04_logs/controller/${run_id}.log"
    mkdir -p "$(dirname "$controller_log")"
    printf 'Controller log: %s\n' "$controller_log"
    "$@" >> "$controller_log" 2>&1 &
    DNA_PIPELINE_ACTIVE_SNAKEMAKE_PID=$!
    progress_interval=60
    [[ -t 1 ]] || progress_interval=300
    local polled_seconds=0
    while kill -0 "$DNA_PIPELINE_ACTIVE_SNAKEMAKE_PID" 2>/dev/null; do
        sleep 2
        polled_seconds=$((polled_seconds + 2))
        kill -0 "$DNA_PIPELINE_ACTIVE_SNAKEMAKE_PID" 2>/dev/null || break
        (( polled_seconds % progress_interval == 0 )) || continue
        last="$(grep -E '^[0-9]+ of [0-9]+ steps' "$controller_log" 2>/dev/null | tail -1)"
        if [[ -n "$last" ]]; then
            done_n="$(printf '%s\n' "$last" | awk '{print $1}')"
            total_n="$(printf '%s\n' "$last" | awk '{print $3}')"
            pct_n="$(printf '%s\n' "$last" | sed -E 's/.*\(([0-9]+)%\).*/\1/')"
            printf '进度: 已完成 %s/%s 个作业（%s%%；完整输出见 %s）\n' \
                "$done_n" "$total_n" "$pct_n" "$controller_log"
        else
            printf '进度: 运行中（完整输出见 %s）\n' "$controller_log"
        fi
    done
    wait "$DNA_PIPELINE_ACTIVE_SNAKEMAKE_PID"
    rc=$?
    DNA_PIPELINE_ACTIVE_SNAKEMAKE_PID=""
    trap - INT TERM HUP
    set -e
    if (( rc != 0 )); then
        echo "ERROR: workflow 执行失败（rc=${rc}）；以下为日志尾部："
        tail -n 40 "$controller_log"
        if grep -q "cannot be locked" "$controller_log" 2>/dev/null; then
            echo "提示: Snakemake 项目锁未释放。确认没有其它 launcher 在本项目上运行后，" >&2
            echo "删除 $DNA_PIPELINE_PROJECT_DIR/.snakemake/locks 再重试。" >&2
        fi
        return "$rc"
    fi
    return 0
}


for internal_helper in \
    "$DNA_PIPELINE_ROOT/core/job_runtime.sh" \
    "$DNA_PIPELINE_ROOT/core/dna_pipeline.py"; do
    [[ -f "$internal_helper" ]] || {
        echo "ERROR: Pipeline helper is missing: $internal_helper" >&2
        exit 2
    }
done

DNA_PIPELINE_WORKFLOW_CONFIG="$DNA_PIPELINE_CONFIG"
DNA_PIPELINE_RUNTIME_TMP=""
if (( DNA_PIPELINE_READ_ONLY == 1 )) || [[ "$DNA_PIPELINE_DATA_MODE" == test ]]; then
    DNA_PIPELINE_RUNTIME_TMP="$(mktemp -d -t dna-pipeline-launch.XXXXXX)"
fi
trap cleanup_launcher_tmp EXIT
trap 'cleanup_launcher_tmp; exit 130' INT
trap 'cleanup_launcher_tmp; exit 143' TERM
trap 'cleanup_launcher_tmp; exit 129' HUP

if [[ "$DNA_PIPELINE_DATA_MODE" == test ]]; then
    [[ "$DNA_PIPELINE_TEST_ROUTE" != auto ]] || DNA_PIPELINE_TEST_ROUTE="hg38-cabernet-bismark"
    DNA_PIPELINE_WORKFLOW_CONFIG="$DNA_PIPELINE_RUNTIME_TMP/effective-config.yaml"
    PYTHONPATH="$DNA_PIPELINE_ROOT/core" "$DNA_PIPELINE_PYTHON" -m dna_pipeline synthesize-config \
        --source "$DNA_PIPELINE_CONFIG" \
        --destination "$DNA_PIPELINE_WORKFLOW_CONFIG" \
        --route "$DNA_PIPELINE_TEST_ROUTE"
fi

DNA_PIPELINE_CONFIG_SCALARS="$(PYTHONPATH="$DNA_PIPELINE_ROOT/core" "$DNA_PIPELINE_PYTHON" \
    -m dna_pipeline config-scalars --config "$DNA_PIPELINE_WORKFLOW_CONFIG")"
eval "$DNA_PIPELINE_CONFIG_SCALARS"
DNA_PIPELINE_METHYLATION_BACKEND="$CFG_METHYLATION_BACKEND"
case "$DNA_PIPELINE_METHYLATION_BACKEND" in
    rastair)
        require_environment_tools "$DNA_PIPELINE_RASTAIR_ENV" rastair "${RASTAIR_REQUIRED[@]}"
        DNA_PIPELINE_RASTAIR_ENV="$(cd "$DNA_PIPELINE_RASTAIR_ENV" && pwd -P)"
        DNA_PIPELINE_RULE_ENV="$DNA_PIPELINE_RASTAIR_ENV"
        ;;
    biscuit)
        require_environment_tools "$DNA_PIPELINE_BISCUIT_ENV" biscuit "${BISCUIT_REQUIRED[@]}"
        DNA_PIPELINE_BISCUIT_ENV="$(cd "$DNA_PIPELINE_BISCUIT_ENV" && pwd -P)"
        DNA_PIPELINE_RULE_ENV="$DNA_PIPELINE_BISCUIT_ENV"
        ;;
    bismark)
        require_environment_tools "$DNA_PIPELINE_BISMARK_ENV" bismark "${BISMARK_RUNTIME_REQUIRED[@]}"
        DNA_PIPELINE_BISMARK_ENV="$(cd "$DNA_PIPELINE_BISMARK_ENV" && pwd -P)"
        DNA_PIPELINE_RULE_ENV="$DNA_PIPELINE_BISMARK_ENV"
        ;;
    *)
        echo "ERROR: analysis.methylation_backend must be biscuit, bismark or rastair." >&2
        exit 2
        ;;
esac
[[ -d "$DNA_PIPELINE_BISCUIT_ENV" ]] && DNA_PIPELINE_BISCUIT_ENV="$(cd "$DNA_PIPELINE_BISCUIT_ENV" && pwd -P)"
[[ -d "$DNA_PIPELINE_BISMARK_ENV" ]] && DNA_PIPELINE_BISMARK_ENV="$(cd "$DNA_PIPELINE_BISMARK_ENV" && pwd -P)"
[[ -d "$DNA_PIPELINE_RASTAIR_ENV" ]] && DNA_PIPELINE_RASTAIR_ENV="$(cd "$DNA_PIPELINE_RASTAIR_ENV" && pwd -P)"
export DNA_PIPELINE_BISCUIT_ENV DNA_PIPELINE_BISMARK_ENV DNA_PIPELINE_RASTAIR_ENV DNA_PIPELINE_RULE_ENV

export PATH="$DNA_PIPELINE_CONTROL_ENV/bin:$PATH"
export CONDA_EXE="$DNA_PIPELINE_CONDA"
unset PYTHONHOME
export PYTHONNOUSERSITE=1
export PYTHONPATH="$DNA_PIPELINE_ROOT/core"

# 一次解析 status JSON 的 status 与 reason 两个字段，避免逐字段 spawn Python。
pipeline_status_fields() {
    "$DNA_PIPELINE_PYTHON" -c 'import json,sys; d=json.loads(sys.argv[1]); print(d.get("status","")); print(d.get("reason",""))' "$1"
}

if (( DNA_PIPELINE_READ_ONLY == 1 )); then
    XDG_CACHE_HOME="$DNA_PIPELINE_RUNTIME_TMP/xdg_cache"
else
    XDG_CACHE_HOME="$DNA_PIPELINE_PROJECT_DIR/05_tmp/xdg_cache"
fi
mkdir -p "$XDG_CACHE_HOME"
export XDG_CACHE_HOME

DNA_PIPELINE_SNAKEMAKE_VERSION="$($DNA_PIPELINE_SNAKEMAKE --version 2>/dev/null | head -n1)"
DNA_PIPELINE_SNAKEMAKE_MAJOR="$(printf '%s' "$DNA_PIPELINE_SNAKEMAKE_VERSION" | sed -E 's/^([0-9]+).*/\1/')"
[[ "$DNA_PIPELINE_SNAKEMAKE_MAJOR" =~ ^[0-9]+$ ]] && (( DNA_PIPELINE_SNAKEMAKE_MAJOR >= 9 )) || {
    echo "ERROR: Snakemake >=9 is required; found $DNA_PIPELINE_SNAKEMAKE_VERSION." >&2
    exit 2
}


# 仅强制无输出的 all 聚合节点，确保仅重发布时也进入持锁成功回调，不强制其科学依赖重算。
dna_pipeline_prepare_common_args() {
    local latency_wait="$1"
    DNA_PIPELINE_COMMON_ARGS=(
        --profile none
        --workflow-profile none
        --keep-going
        --retries 1
        --rerun-incomplete
        --forcerun all
        --rerun-triggers mtime input params code software-env
        --latency-wait "$latency_wait"
        --show-failed-logs
        --software-deployment-method conda
        --conda-frontend conda
    )
}

dna_pipeline_prepare_run_snapshot() {
    local -a snapshot_args
    snapshot_args=(
        --pipeline-root "$DNA_PIPELINE_ROOT"
        --project "$DNA_PIPELINE_PROJECT_DIR"
        --config "$DNA_PIPELINE_WORKFLOW_CONFIG"
        --control-env "$DNA_PIPELINE_CONTROL_ENV"
        --biscuit-env "$DNA_PIPELINE_BISCUIT_ENV"
        --bismark-env "$DNA_PIPELINE_BISMARK_ENV"
        --rastair-env "$DNA_PIPELINE_RASTAIR_ENV"
    )
    if [[ "$DNA_PIPELINE_DATA_MODE" == test ]]; then
        snapshot_args+=(--allow-external-config)
    fi
    DNA_PIPELINE_RUN_SNAPSHOT="$($DNA_PIPELINE_PYTHON -m dna_pipeline snapshot-create "${snapshot_args[@]}")"
    [[ -f "$DNA_PIPELINE_RUN_SNAPSHOT" ]] || {
        echo "ERROR: run-start snapshot was not created: $DNA_PIPELINE_RUN_SNAPSHOT" >&2
        exit 2
    }
    SNAPSHOT_CONFIG="$(dirname "$DNA_PIPELINE_RUN_SNAPSHOT")/inputs/config.yaml"
    [[ -f "$SNAPSHOT_CONFIG" ]] || {
        echo "ERROR: run snapshot has no config copy: $SNAPSHOT_CONFIG" >&2
        exit 2
    }
    DNA_PIPELINE_DELIVERY_ID="$(basename "$(dirname "$DNA_PIPELINE_RUN_SNAPSHOT")")"
    export DNA_PIPELINE_DELIVERY_ID
    export DNA_PIPELINE_RUN_SNAPSHOT
}

dna_pipeline_prepare_read_only_inputs() {
    SNAPSHOT_CONFIG="$DNA_PIPELINE_WORKFLOW_CONFIG"
    DNA_PIPELINE_DELIVERY_ID="dryrun-$$"
    export DNA_PIPELINE_DELIVERY_ID
}

if (( DNA_PIPELINE_READ_ONLY == 0 )); then
    dna_pipeline_prepare_run_snapshot
    cleanup_launcher_tmp
    trap cleanup_launcher_tmp EXIT
fi
DNA_PIPELINE_STATUS_ARGS=(--project "$DNA_PIPELINE_PROJECT_DIR")
if [[ -n "${DNA_PIPELINE_RUN_SNAPSHOT:-}" ]]; then
    DNA_PIPELINE_STATUS_ARGS+=(--current-snapshot "$DNA_PIPELINE_RUN_SNAPSHOT")
fi
DNA_PIPELINE_STATUS_JSON="$($DNA_PIPELINE_PYTHON -m dna_pipeline project-status "${DNA_PIPELINE_STATUS_ARGS[@]}")"
DNA_PIPELINE_STATUS_LINES="$(pipeline_status_fields "$DNA_PIPELINE_STATUS_JSON")"
DNA_PIPELINE_PROJECT_STATUS="$(printf '%s\n' "$DNA_PIPELINE_STATUS_LINES" | sed -n 1p)"
DNA_PIPELINE_PROJECT_STATUS_REASON="$(printf '%s\n' "$DNA_PIPELINE_STATUS_LINES" | sed -n 2p)"

if [[ "$DNA_PIPELINE_PROJECT_STATUS" == complete_current ]]; then
    printf 'Alopex delivery is complete and current: %s\n' \
        "$DNA_PIPELINE_PROJECT_DIR/03_results/run_manifest.json"
    printf 'No workflow work is required. Snakemake state was left untouched.\n'
    "$DNA_PIPELINE_PYTHON" -m dna_pipeline clean-demux-scratch --project "$DNA_PIPELINE_PROJECT_DIR"
    exit 0
elif [[ "$DNA_PIPELINE_PROJECT_STATUS" == resumable ]]; then
    if [[ -n "$DNA_PIPELINE_PROJECT_STATUS_REASON" ]]; then
        printf 'Resuming existing project state: %s\n' "$DNA_PIPELINE_PROJECT_STATUS_REASON"
    else
        printf 'Resuming existing project state; Snakemake will decide the minimal rerun set.\n'
    fi
fi

# 按平台选择执行器，已有 Slurm allocation 仍提交规则作业。
detect_executor() {
    local system_name
    system_name="$(uname -s)"

    case "$system_name" in
        Darwin) printf '%s\n' local ;;
        Linux)
            if [[ -n "$DNA_PIPELINE_DOCTOR_LOCAL_CORES" ]]; then
                printf '%s\n' local
            else
                command -v sbatch >/dev/null 2>&1 || { echo "ERROR: Linux execution requires sbatch." >&2; return 127; }
                command -v squeue >/dev/null 2>&1 || { echo "ERROR: Linux execution requires squeue." >&2; return 127; }
                printf '%s\n' slurm
            fi ;;
        *) echo "ERROR: unsupported OS: $system_name" >&2; return 2 ;;
    esac
}
DNA_PIPELINE_EXECUTOR="$(detect_executor)"
DNA_PIPELINE_IN_SLURM_ALLOCATION=0
if [[ "$(uname -s)" == Linux && -n "${SLURM_JOB_ID:-}" ]]; then
    DNA_PIPELINE_IN_SLURM_ALLOCATION=1
fi
export DNA_PIPELINE_EXECUTOR DNA_PIPELINE_IN_SLURM_ALLOCATION DNA_PIPELINE_DATA_MODE DNA_PIPELINE_TEST_ROUTE

# 本地执行器的 mem_mb 调度预算：SLURM allocation 内取作业实际内存，否则取宿主总内存；
# 使多个整基因组 index 级作业按声明内存串行，避免并发驻留击穿 cgroup/物理内存。
dna_pipeline_local_mem_budget_mib() {
    local budget=""
    if [[ "${SLURM_MEM_PER_NODE:-}" =~ ^[1-9][0-9]*$ ]]; then
        budget="$SLURM_MEM_PER_NODE"
    elif [[ "${SLURM_MEM_PER_CPU:-}" =~ ^[1-9][0-9]*$ && "${SLURM_CPUS_ON_NODE:-}" =~ ^[1-9][0-9]*$ ]]; then
        budget=$(( SLURM_MEM_PER_CPU * SLURM_CPUS_ON_NODE ))
    elif [[ "$(uname -s)" == Darwin ]]; then
        budget=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1048576 ))
    else
        budget=$(( $(awk '/^MemTotal:/ {print $2}' /proc/meminfo 2>/dev/null || echo 0) / 1024 ))
    fi
    [[ "$budget" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: 无法确定本地 mem_mb 预算。" >&2; return 2; }
    echo "$budget"
}

# 在 macOS 或 Doctor 内部小型路线运行本地执行器。
run_local() {
    local cores latency_wait protocol reference execution_label mem_budget
    cores="$CFG_LOCAL_CORES"
    [[ -z "$DNA_PIPELINE_DOCTOR_LOCAL_CORES" ]] || cores="$DNA_PIPELINE_DOCTOR_LOCAL_CORES"
    execution_label="local"
    latency_wait="$CFG_LATENCY_WAIT"
    protocol="$CFG_PROTOCOL"
    reference="$CFG_SPECIES"
    [[ "$cores" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: runtime.local_executor_cores must be positive." >&2; return 2; }
    [[ "$latency_wait" =~ ^[0-9]+$ ]] || { echo "ERROR: runtime.latency_wait_seconds must be >=0." >&2; return 2; }
    mem_budget="$(dna_pipeline_local_mem_budget_mib)" || return $?
    dna_pipeline_prepare_common_args "$latency_wait"
    if (( DNA_PIPELINE_READ_ONLY == 1 )); then
        dna_pipeline_prepare_read_only_inputs
    fi
    echo "Starting Alopex with local executor..."
    echo "Config: $DNA_PIPELINE_WORKFLOW_CONFIG"
    echo "Project: $DNA_PIPELINE_PROJECT_DIR"
    echo "Executor: $execution_label"
    echo "Protocol: $protocol"
    echo "Backend: $DNA_PIPELINE_METHYLATION_BACKEND"
    echo "Reference: $reference"
    echo "Cores: $cores"
    echo "Mem budget: ${mem_budget} MiB"
    export DNA_PIPELINE_LOCAL_CORES="$cores"
    if (( DNA_PIPELINE_READ_ONLY == 0 )); then echo "Run snapshot: $DNA_PIPELINE_RUN_SNAPSHOT"; fi
    if (( DNA_PIPELINE_READ_ONLY == 1 )); then
        cd "$DNA_PIPELINE_RUNTIME_TMP"
    else
        cd "$DNA_PIPELINE_PROJECT_DIR"
    fi
    DNA_PIPELINE_CMD=("$DNA_PIPELINE_SNAKEMAKE" --snakefile "$DNA_PIPELINE_ROOT/core/Snakefile" --executor local --cores "$cores" --resources mem_mb="$mem_budget" "${DNA_PIPELINE_COMMON_ARGS[@]}")
    (( DNA_PIPELINE_READ_ONLY == 0 )) || DNA_PIPELINE_CMD+=(--dry-run)
    DNA_PIPELINE_CMD+=(--configfile "$SNAPSHOT_CONFIG")
    local controller_run_id
    controller_run_id="$DNA_PIPELINE_DELIVERY_ID"
    run_snakemake_tracked "$controller_run_id" "${DNA_PIPELINE_CMD[@]}"
}

# 在当前 shell 运行 Slurm controller，规则作业通过 sbatch 提交。
run_slurm() {
    local jobs controller_cores latency_wait account partition qos configured_tmp project_hash plugin_version protocol effective_account
    local -a default_resources plugin_args cmd
    command -v sbatch >/dev/null 2>&1 || { echo "ERROR: sbatch not found." >&2; return 127; }
    command -v squeue >/dev/null 2>&1 || { echo "ERROR: squeue not found." >&2; return 127; }
    local logdir
    if (( DNA_PIPELINE_READ_ONLY == 1 )); then
        logdir="$DNA_PIPELINE_RUNTIME_TMP/slurm"
    else
        logdir="$DNA_PIPELINE_PROJECT_DIR/04_logs/slurm"
    fi
    mkdir -p "$logdir"
    jobs="$CFG_SLURM_JOBS"; controller_cores="$CFG_SLURM_CONTROLLER_CORES"; latency_wait="$CFG_LATENCY_WAIT"
    account="$CFG_SLURM_ACCOUNT"; partition="$CFG_SLURM_PARTITION"; qos="$CFG_SLURM_QOS"; configured_tmp="$CFG_SLURM_NODE_TMPDIR"
    [[ "$jobs" =~ ^[1-9][0-9]*$ && "$controller_cores" =~ ^[1-9][0-9]*$ && "$latency_wait" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid Slurm/runtime numeric settings." >&2; return 2; }
    if [[ -n "$configured_tmp" ]]; then
        DNA_PIPELINE_NODE_TMPDIR="$(PYTHONPATH="$DNA_PIPELINE_ROOT/core" "$DNA_PIPELINE_PYTHON" \
            -m dna_pipeline resolve-tmpdir --path "$configured_tmp")"
        export DNA_PIPELINE_NODE_TMPDIR
    else
        unset DNA_PIPELINE_NODE_TMPDIR
    fi
    project_hash="$(project_hash8)"
    plugin_version="$($DNA_PIPELINE_PYTHON -c 'from importlib.metadata import version;print(version("snakemake-executor-plugin-slurm"))' 2>/dev/null || true)"
    [[ -n "$plugin_version" ]] || { echo "ERROR: SLURM executor plugin is missing." >&2; return 127; }
    dna_pipeline_prepare_common_args "$latency_wait"
    default_resources=(); effective_account="$account"
    if [[ -n "$partition" ]]; then default_resources+=("slurm_partition=$partition"); fi
    if [[ -n "$effective_account" ]]; then default_resources+=("slurm_account=$effective_account"); fi
    plugin_args=(--slurm-jobname-prefix "Alopex_${project_hash}" --slurm-logdir "$logdir" --slurm-delete-logfiles-older-than 0 --slurm-status-command squeue)
    if [[ -z "$effective_account" ]]; then plugin_args+=(--slurm-no-account); fi
    if [[ -n "$qos" ]]; then plugin_args+=(--slurm-qos "$qos"); fi
    if (( DNA_PIPELINE_READ_ONLY == 1 )); then
        dna_pipeline_prepare_read_only_inputs
    fi
    protocol="$CFG_PROTOCOL"
    reference="$CFG_SPECIES"
    echo "Starting Alopex with Snakemake SLURM executor..."
    echo "Config: $DNA_PIPELINE_WORKFLOW_CONFIG"; echo "Project: $DNA_PIPELINE_PROJECT_DIR"
    if (( DNA_PIPELINE_IN_SLURM_ALLOCATION == 1 )); then
        echo "Executor: slurm (plugin $plugin_version; controller=current allocation ${SLURM_JOB_ID})"
    else
        echo "Executor: slurm (plugin $plugin_version; controller=current login shell)"
    fi
    echo "Protocol: $protocol"; echo "Backend: $DNA_PIPELINE_METHYLATION_BACKEND"; echo "Reference: $reference"; echo "SLURM jobs: $jobs (controller cores=$controller_cores)"
    if (( DNA_PIPELINE_READ_ONLY == 0 )); then echo "Run snapshot: $DNA_PIPELINE_RUN_SNAPSHOT"; fi
    if (( DNA_PIPELINE_READ_ONLY == 1 )); then
        cd "$DNA_PIPELINE_RUNTIME_TMP"
    else
        cd "$DNA_PIPELINE_PROJECT_DIR"
    fi
    cmd=("$DNA_PIPELINE_SNAKEMAKE" --snakefile "$DNA_PIPELINE_ROOT/core/Snakefile" --executor slurm --max-jobs-per-timespan 30/1s --max-status-checks-per-second 5 --scheduler greedy --jobs "$jobs" --local-cores "$controller_cores" "${DNA_PIPELINE_COMMON_ARGS[@]}")
    if (( ${#default_resources[@]} > 0 )); then cmd+=(--default-resources "${default_resources[@]}"); fi
    cmd+=("${plugin_args[@]}")
    (( DNA_PIPELINE_READ_ONLY == 0 )) || cmd+=(--dry-run)
    cmd+=(--configfile "$SNAPSHOT_CONFIG")
    local controller_run_id
    controller_run_id="$DNA_PIPELINE_DELIVERY_ID"
    run_snakemake_tracked "$controller_run_id" "${cmd[@]}"
}

case "$DNA_PIPELINE_EXECUTOR" in
    local) run_local ;;
    slurm) run_slurm ;;
esac

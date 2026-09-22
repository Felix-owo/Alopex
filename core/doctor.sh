#!/usr/bin/env bash

DOCTOR_SELF_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/$(basename "${BASH_SOURCE[0]}")"
PIPELINE_ROOT="$(cd "$(dirname "$DOCTOR_SELF_PATH")/.." && pwd -P)"
if [[ ! -f "$PIPELINE_ROOT/core/doctor.sh" && -n "${DNA_PIPELINE_ROOT:-}" && -f "$DNA_PIPELINE_ROOT/core/doctor.sh" ]]; then
    PIPELINE_ROOT="$DNA_PIPELINE_ROOT"
fi
export DNA_PIPELINE_ROOT="$PIPELINE_ROOT"


TXN_LOCK_HELD=0

# 从两列 TSV（key\tvalue）读取首个匹配键的值；文件缺失返回非零。
txn_read_tsv_field() {
    local file="$1" field="$2"
    [[ -f "$file" ]] || return 1
    awk -F '\t' -v key="$field" '$1 == key {print $2; exit}' "$file"
}

txn_lock_hostname() {
    hostname 2>/dev/null || uname -n
}

txn_read_owner_field() {
    txn_read_tsv_field "$1/owner.tsv" "$2"
}

# 用 squeue 判断锁 owner 是否存活；作业终结或号码复用才判死，证据不足保持 unknown。
txn_slurm_owner_alive() {
    local job_id="$1" lock_started_epoch="$2"
    local squeue_output row state start_epoch
    command -v squeue >/dev/null 2>&1 || { printf 'unknown'; return 0; }
    if ! squeue_output="$(SLURM_TIME_FORMAT='%s' squeue -h -j "$job_id" -o '%T|%S' 2>&1)"; then
        if [[ "$squeue_output" == *"Invalid job id"* ]]; then
            printf 'dead'
        else
            printf 'unknown'
        fi
        return 0
    fi
    row="$(head -n 1 <<<"$squeue_output")"
    [[ -n "$row" ]] || { printf 'dead'; return 0; }
    state="${row%%|*}"
    start_epoch="${row#*|}"
    case "$state" in
        COMPLETED|FAILED|CANCELLED*|TIMEOUT|OUT_OF_MEMORY|BOOT_FAIL|NODE_FAIL|PREEMPTED*|DEADLINE|REVOKED)
            printf 'dead'
            ;;
        PENDING*|CONFIGURING)
            printf 'unknown'
            ;;
        RUNNING)
            [[ "$start_epoch" =~ ^[0-9]+$ && "$lock_started_epoch" =~ ^[0-9]+$ ]] || { printf 'unknown'; return 0; }
            if (( start_epoch <= lock_started_epoch + 60 )); then
                printf 'alive'
            else
                printf 'dead'
            fi
            ;;
        *)
            printf 'unknown'
            ;;
    esac
    return 0
}

# 仅按同主机死 PID 或 Slurm 作业死亡证据回收锁；无法证明时保持 busy。
txn_try_remove_stale_lock() {
    local lock_dir="$1" owner_host owner_pid owner_job owner_started_epoch alive
    [[ -d "$lock_dir" ]] || return 1
    owner_host="$(txn_read_owner_field "$lock_dir" host 2>/dev/null || true)"
    owner_pid="$(txn_read_owner_field "$lock_dir" pid 2>/dev/null || true)"
    [[ -n "$owner_host" && "$owner_pid" =~ ^[0-9]+$ ]] || return 1
    if [[ "$owner_host" == "$(txn_lock_hostname)" ]]; then
        kill -0 "$owner_pid" 2>/dev/null && return 1
        echo "WARNING: 回收已失效的事务锁： $lock_dir (host=$owner_host pid=$owner_pid)" >&2
        rm -rf -- "$lock_dir"
        return 0
    fi
    owner_job="$(txn_read_owner_field "$lock_dir" slurm_job 2>/dev/null || true)"
    [[ "$owner_job" =~ ^[0-9]+$ ]] || return 1
    owner_started_epoch="$(txn_read_owner_field "$lock_dir" started_epoch 2>/dev/null || true)"
    alive="$(txn_slurm_owner_alive "$owner_job" "$owner_started_epoch")"
    [[ "$alive" == "dead" ]] || return 1
    echo "WARNING: 回收已失效的事务锁： $lock_dir (host=$owner_host pid=$owner_pid slurm_job=$owner_job)" >&2
    rm -rf -- "$lock_dir"
    return 0
}

# mkdir 原子占位；成功后写 owner.tsv 并置 TXN_LOCK_HELD=1；忙锁返回 75。
txn_acquire_lock() {
    local lock_dir="$1" label="${2:-transaction}" attempt owner_host owner_pid owner_job job_state_line
    mkdir -p "$(dirname "$lock_dir")"
    for attempt in 1 2; do
        if mkdir "$lock_dir" 2>/dev/null; then
            {
                printf 'pid\t%s\n' "$$"
                printf 'host\t%s\n' "$(txn_lock_hostname)"
                printf 'started_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
                printf 'started_epoch\t%s\n' "$(date +%s)"
                [[ "${SLURM_JOB_ID:-}" =~ ^[0-9]+$ ]] && printf 'slurm_job\t%s\n' "$SLURM_JOB_ID"
                printf 'command\t%s\n' "$0"
            } > "$lock_dir/owner.tsv"
            TXN_LOCK_HELD=1
            return 0
        fi
        if (( attempt == 1 )) && txn_try_remove_stale_lock "$lock_dir"; then
            continue
        fi
        owner_host="$(txn_read_owner_field "$lock_dir" host 2>/dev/null || printf '?')"
        owner_pid="$(txn_read_owner_field "$lock_dir" pid 2>/dev/null || printf '?')"
        owner_job="$(txn_read_owner_field "$lock_dir" slurm_job 2>/dev/null || true)"
        echo "ERROR: 已有另一项 $label 正在运行。" >&2
        echo "  lock : $lock_dir" >&2
        echo "  owner: host=$owner_host pid=$owner_pid${owner_job:+ slurm_job=$owner_job}" >&2
        if [[ "$owner_job" =~ ^[0-9]+$ ]] && command -v squeue >/dev/null 2>&1; then
            job_state_line="$(squeue -h -j "$owner_job" -o '%T on %N (%M)' 2>/dev/null | head -n 1)"
            [[ -n "$job_state_line" ]] && echo "  job  : $owner_job -> $job_state_line" >&2
        fi
        echo "  hint : 若确认没有进行中的构建/修复，可删除该锁目录后重试" >&2
        return 75
    done
    return 75
}

# 只有「同主机 + PID == $$」的 owner 才允许删除锁目录。
txn_release_lock() {
    local lock_dir="$1" owner_host owner_pid
    (( TXN_LOCK_HELD == 1 )) || return 0
    [[ -d "$lock_dir" ]] || { TXN_LOCK_HELD=0; return 0; }
    owner_host="$(txn_read_owner_field "$lock_dir" host 2>/dev/null || true)"
    owner_pid="$(txn_read_owner_field "$lock_dir" pid 2>/dev/null || true)"
    if [[ "$owner_host" == "$(txn_lock_hostname)" && "$owner_pid" == "$$" ]]; then
        rm -rf -- "$lock_dir"
    else
        echo "WARNING: 当前进程不是锁 owner，拒绝删除： $lock_dir" >&2
    fi
    TXN_LOCK_HELD=0
}

txn_pid_is_running() {
    kill -0 "$1" 2>/dev/null
}

# 轮询等待 PID 退出，超过宽限秒数返回非零；进程为本 shell 子进程时一并 reap。
txn_wait_for_pid_exit() {
    local pid="$1" grace_seconds="$2" deadline now
    deadline=$(($(date +%s)+grace_seconds))
    while txn_pid_is_running "$pid"; do
        now="$(date +%s)"
        (( now < deadline )) || return 1
        sleep 0.1
    done
    wait "$pid" 2>/dev/null || true
    return 0
}

# 对进程树发信号；KILL 时先递归展开子进程再杀父进程。
txn_signal_process_tree() {
    local pid="$1" signal_name="$2" child
    if [[ "$signal_name" == KILL ]] && command -v pgrep >/dev/null 2>&1; then
        for child in $(pgrep -P "$pid" 2>/dev/null || true); do
            [[ "$child" =~ ^[0-9]+$ ]] || continue
            txn_signal_process_tree "$child" "$signal_name"
        done
    fi
    kill -s "$signal_name" "$pid" 2>/dev/null || true
}

# TERM -> 宽限等待 -> KILL 的有界终止；日志走 stderr，调用方可整体重定向。
txn_terminate_pid_bounded() {
    local pid="$1" label="$2" grace_seconds="${3:-15}" term_already_sent="${4:-0}"
    [[ "$grace_seconds" =~ ^[1-9][0-9]*$ ]] || grace_seconds=15
    [[ "$pid" =~ ^[0-9]+$ ]] || return 0
    if ! txn_pid_is_running "$pid"; then
        wait "$pid" 2>/dev/null || true
        return 0
    fi
    if (( term_already_sent == 0 )); then
        printf '向 %s（pid=%s）发送 TERM\n' "$label" "$pid" >&2
        txn_signal_process_tree "$pid" TERM
    fi
    if txn_wait_for_pid_exit "$pid" "$grace_seconds"; then
        return 0
    fi
    printf '%s 未在 %ss 内退出，发送 KILL（pid=%s）\n' \
        "$label" "$grace_seconds" "$pid" >&2
    txn_signal_process_tree "$pid" KILL
    if ! txn_wait_for_pid_exit "$pid" 2; then
        printf 'WARNING: %s pid=%s 在 KILL 后仍未退出。\n' "$label" "$pid" >&2
    fi
}

CONTROL_RUNTIME_REQUIRED=(
    bin/python
    bin/snakemake
    bin/conda
)

CONTROL_REQUIRED=(
    "${CONTROL_RUNTIME_REQUIRED[@]}"
    bin/cargo
    bin/rustc
)

BISCUIT_REQUIRED=(
    bin/python
    bin/biscuit
    bin/dupsifter
    bin/cutadapt
    bin/samtools
    bin/bgzip
    bin/tabix
    bin/awk
    bin/multiqc
    bin/zstd
    bin/gzip
)

BISMARK_RUNTIME_REQUIRED=(
    bin/umi_tools
    bin/python
    bin/bismark
    bin/bowtie2
    bin/cutadapt
    bin/samtools
    bin/awk
    bin/multiqc
    bin/zstd
    bin/gzip
)

BISMARK_REQUIRED=(
    "${BISMARK_RUNTIME_REQUIRED[@]}"
    bin/bismark_genome_preparation
    bin/bowtie2-build
)

RASTAIR_REQUIRED=(bin/python bin/bwa bin/rastair bin/cutadapt bin/samtools bin/awk bin/multiqc bin/zstd bin/gzip)

NOTEBOOK_REQUIRED=(
    bin/python
    bin/jupyter
    bin/jupyter-lab
)

doctor_usage() {
    cat <<'EOF'
Alopex Doctor

用法：
  bash core/doctor.sh                      检查/修复并执行三路线回归
  bash core/doctor.sh <conda|reference|demux> [--check]   单段构建/检查
  bash core/doctor.sh download             只预热下载缓存（conda 包 + cargo 依赖，需网络）
EOF
}





CONDA_ROOT="$PIPELINE_ROOT/conda"

RELEASE_STORE="$CONDA_ROOT/releases"

CURRENT_RELEASE="$CONDA_ROOT/current"

PACKAGE_CACHE="$CONDA_ROOT/package_cache"

SETUP_HOME="$CONDA_ROOT/setup_home"
SETUP_ENVS="$CONDA_ROOT/setup_envs"
SETUP_CACHE="$CONDA_ROOT/setup_cache"
SETUP_CONFIG="$CONDA_ROOT/setup_config"

CLEAN_ENVIRONMENT=(
    -u CONDA_PREFIX
    -u CONDA_PREFIX_1
    -u CONDA_DEFAULT_ENV
    -u CONDA_PROMPT_MODIFIER
    -u CONDA_SHLVL
    -u CONDA_EXE
    -u CONDA_PYTHON_EXE
    -u CONDA_ROOT_PREFIX
    -u CONDARC

    -u CONDA_SUBDIR
    -u CONDA_OVERRIDE_ARCHSPEC
    -u CONDA_OVERRIDE_CUDA
    -u CONDA_OVERRIDE_GLIBC
    -u CONDA_OVERRIDE_LINUX
    -u _CE_CONDA
    -u _CE_M

    -u PYTHONHOME
    -u PYTHONPATH
    -u PYTHONSTARTUP
    -u PYTHONUSERBASE
    -u VIRTUAL_ENV

    -u LD_LIBRARY_PATH
    -u LD_PRELOAD
    -u LIBRARY_PATH
    -u PKG_CONFIG_PATH
    -u DYLD_LIBRARY_PATH
    -u CMAKE_PREFIX_PATH
    -u CONDA_BUILD_SYSROOT
    -u CONDA_TOOLCHAIN_BUILD
    -u CONDA_TOOLCHAIN_HOST
    -u JAVA_LD_LIBRARY_PATH

    HOME="$SETUP_HOME"
    XDG_CACHE_HOME="$SETUP_CACHE"
    XDG_CONFIG_HOME="$SETUP_CONFIG"

    CONDA_PKGS_DIRS="$PACKAGE_CACHE"
    CONDA_ENVS_PATH="$SETUP_ENVS"
)

ENVIRONMENT_LOCK="$CONDA_ROOT/.locks/environments.lock.d"

ENVS_SPEC="$PIPELINE_ROOT/core/envs.yaml"
ENVS_BUILD_DIR="$CONDA_ROOT/.build/env-specs"

usage_conda() {
    cat <<'USAGE'
用法：
  bash core/doctor.sh conda [--check]

作用：
  自动寻找系统中已经存在、能够工作的 Conda 作为 Seed Conda；
  Seed Conda 只负责创建 control candidate。
  control 创建后必须包含 Conda >=26，随后 biscuit（含 Dupsifter）/ bismark / rastair / notebook 均由 control Conda 创建。
  变化的环境才重新创建；未变化且可用的角色通过 release 内 symlink 复用。全部验证通过后，再原子切换 conda/current。

Seed Conda 自动搜索顺序：
  1. 当前 Pipeline current/control 中的 conda；
  2. 当前 shell 的 $CONDA_PREFIX/bin/conda；
  3. PATH 中的 conda executable。

如果这些位置都找不到可用 Conda：直接报错退出。
本脚本不会下载、安装或 bootstrap Miniforge/Conda。
也不接受手工 Conda 路径覆盖参数。

参数：
  --check
      只读验证 conda/current、release manifest、5 个 environment identity、
      YAML spec SHA、native platform 与五套环境 smoke test；不构建、不修复。

  -h, --help
      显示本帮助。

环境状态判断、自动修复和其它系统检查请使用：
  bash core/doctor.sh
USAGE
}

CHECK_ONLY=0

HOST_PATH="/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"




run_in_sanitized_environment() {
    local executable_dir="$1"
    shift
    env \
        "${CLEAN_ENVIRONMENT[@]}" \
        PATH="$executable_dir:$HOST_PATH" \
        PYTHONNOUSERSITE=1 \
        "$@"
}

native_conda_subdir() {
    local system machine
    system="$(uname -s)"
    machine="$(uname -m)"

    case "$system/$machine" in
        Linux/x86_64|Linux/amd64)
            printf 'linux-64\n'
            ;;
        Linux/aarch64|Linux/arm64)
            printf 'linux-aarch64\n'
            ;;
        Darwin/arm64)
            printf 'osx-arm64\n'
            ;;
        Darwin/x86_64)
            printf 'osx-64\n'
            ;;
        *)
            echo "ERROR: 当前平台暂不支持：$system/$machine" >&2
            return 2
            ;;
    esac
}
NATIVE_SUBDIR="$(native_conda_subdir)"


# env-spec 使用本次 control、现有 control 或已验证 Seed 的 Python，不依赖宿主 PATH。
envs_spec_python() {
    if [[ -x "${CONDA_PYTHON_BIN:-}" ]]; then
        printf '%s\n' "$CONDA_PYTHON_BIN"
    elif [[ -x "$CONDA_ROOT/current/control/bin/python" ]]; then
        printf '%s\n' "$CONDA_ROOT/current/control/bin/python"
    else
        seed_python_path "$SEED_CONDA_DIR"
    fi
}

# 抽取角色 section 为独立 conda spec 文件（conda env create 的输入）。
extract_env_spec() {
    local role="$1" dst python_bin
    dst="$ENVS_BUILD_DIR/$role.yaml"
    mkdir -p "$ENVS_BUILD_DIR" || return $?
    python_bin="$(envs_spec_python)" || return $?
    PYTHONPATH="$PIPELINE_ROOT/core" "$python_bin" -m dna_pipeline env-spec \
        --envs "$ENVS_SPEC" --role "$role" --destination "$dst" >&2 || return $?
    printf '%s\n' "$dst"
}

require_specification() {
    [[ -f "$1" ]] || {
        echo "ERROR: 缺少 Conda 环境定义文件：$1" >&2
        return 2
    }
}


prefix_has_tools() {
    local prefix="$1" relative
    shift

    for relative in "$@"; do
        [[ -x "$prefix/$relative" ]] || return 1
    done
}

prefix_python_imports() {
    local prefix="$1"
    shift

    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/python" - "$@" <<'PY' >/dev/null 2>&1
import importlib
import sys

for module_name in sys.argv[1:]:
    importlib.import_module(module_name)
PY
}

# noarch 包不影响判断；只要发现 linux-64/osx-arm64 等非 noarch 包来自其他平台就失败。
prefix_matches_native_platform() {
    local prefix="$1"
    local python="$prefix/bin/python"

    [[ -x "$python" ]] || return 1

    run_in_sanitized_environment "$prefix/bin" "$python" - "$prefix" "$NATIVE_SUBDIR" <<'PY' >/dev/null 2>&1
import json
import pathlib
import sys

prefix = pathlib.Path(sys.argv[1])
expected = sys.argv[2]
meta_dir = prefix / "conda-meta"
seen_native = False

for path in meta_dir.glob("*.json"):
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        raise SystemExit(1)
    subdir = payload.get("subdir")
    if not subdir or subdir == "noarch":
        continue
    seen_native = True
    if subdir != expected:
        raise SystemExit(1)

raise SystemExit(0 if seen_native else 1)
PY
}

# 如果系统存在 GNU timeout，则最多等待 30 秒；macOS 默认没有 timeout 时直接运行。
run_conda_probe() {
    local executable="$1" executable_dir
    shift
    executable_dir="$(cd "$(dirname "$executable")" && pwd -P)"

    if command -v timeout >/dev/null 2>&1; then
        run_in_sanitized_environment "$executable_dir" timeout 30s "$executable" "$@"
    else
        run_in_sanitized_environment "$executable_dir" "$executable" "$@"
    fi
}

seed_python_path() {
    local conda_dir="$1" python
    python="$conda_dir/python"

    if [[ ! -x "$python" && "$(basename "$conda_dir")" == "condabin" ]]; then
        python="$(cd "$conda_dir/../bin" && pwd -P)/python"
    fi

    [[ -x "$python" ]] || return 1
    printf '%s\n' "$python"
}

validate_seed_conda_path() {
    local candidate="$1" version_output

    [[ "$candidate" == /* ]] || return 1
    [[ "$(basename "$candidate")" == "conda" ]] || return 1
    [[ -f "$candidate" && -x "$candidate" ]] || return 1

    version_output="$(run_conda_probe "$candidate" --version 2>&1)" || return 1
    printf '%s\n' "$version_output" | grep -Eq '^conda[[:space:]]+[0-9]+' || return 1

    run_conda_probe "$candidate" info --json >/dev/null 2>&1 || return 1

    printf '%s/%s\n' \
        "$(cd "$(dirname "$candidate")" && pwd -P)" \
        "$(basename "$candidate")"
}
resolve_seed_conda() {
    local candidate="" resolved=""

    candidate="$CURRENT_RELEASE/control/bin/conda"
    resolved="$(validate_seed_conda_path "$candidate" 2>/dev/null || true)"
    if [[ -n "$resolved" ]]; then
        echo "使用 Pipeline 当前 control Conda 作为 Seed：$resolved" >&2
        printf '%s\n' "$resolved"
        return 0
    fi

    if [[ -n "${CONDA_PREFIX:-}" && -e "$CONDA_PREFIX/bin/conda" ]]; then
        candidate="$CONDA_PREFIX/bin/conda"
        resolved="$(validate_seed_conda_path "$candidate" 2>/dev/null || true)"
        if [[ -n "$resolved" ]]; then
            echo "使用当前 CONDA_PREFIX 中的 Conda 作为 Seed：$resolved" >&2
            printf '%s\n' "$resolved"
            return 0
        fi
    fi

    candidate="$(type -P conda 2>/dev/null || true)"
    if [[ -n "$candidate" ]]; then
        candidate="$(cd "$(dirname "$candidate")" 2>/dev/null && printf '%s/%s\n' "$(pwd -P)" "$(basename "$candidate")")" || candidate=""
        if [[ -n "$candidate" ]]; then
            resolved="$(validate_seed_conda_path "$candidate" 2>/dev/null || true)"
            if [[ -n "$resolved" ]]; then
                echo "使用 PATH 中的 Conda 作为 Seed：$resolved" >&2
                printf '%s\n' "$resolved"
                return 0
            fi
        fi
    fi

    echo "ERROR: 找不到任何可工作的 Conda executable。" >&2
    echo "  core/doctor.sh conda 已自动检查：" >&2
    echo "    1. conda/current/control/bin/conda" >&2
    echo "    2. \$CONDA_PREFIX/bin/conda（若当前 shell 已激活 Conda）" >&2
    echo "    3. PATH 中的 conda executable" >&2
    echo "  本脚本不安装 Conda。请先自行安装/配置一个可工作的 Conda，再重新运行。" >&2
    return 127
}

validate_seed_context() {
    local info_json seed_python

    seed_python="$(seed_python_path "$SEED_CONDA_DIR")" || {
        echo "ERROR: Seed Conda 没有 prefix-local Python。" >&2
        return 1
    }

    if command -v timeout >/dev/null 2>&1; then
        info_json="$(run_in_sanitized_environment "$SEED_CONDA_DIR" timeout 120s \
            "$SEED_CONDA" info --json)" || {
            echo "ERROR: 无法读取 Seed Conda context。" >&2
            return 1
        }
    else
        info_json="$(run_in_sanitized_environment "$SEED_CONDA_DIR" \
            "$SEED_CONDA" info --json)" || {
            echo "ERROR: 无法读取 Seed Conda context。" >&2
            return 1
        }
    fi

    printf '%s\n' "$info_json" | run_in_sanitized_environment "$SEED_CONDA_DIR" \
        "$seed_python" -c '
import json
import sys

expected_cache = sys.argv[1]
expected_subdir = sys.argv[2]
payload = json.load(sys.stdin)

actual_cache = payload.get("pkgs_dirs")
if actual_cache != [expected_cache]:
    raise SystemExit(
        "ERROR: Seed Conda package cache 超出 Pipeline 私有目录："
        f"expected={[expected_cache]!r}, actual={actual_cache!r}"
    )

actual_subdir = payload.get("subdir")
if actual_subdir and actual_subdir != expected_subdir:
    raise SystemExit(
        "ERROR: Seed Conda 平台不匹配："
        f"expected={expected_subdir!r}, actual={actual_subdir!r}"
    )

print(f"Seed package cache: {expected_cache}", file=sys.stderr)
print(f"Native Conda subdir: {expected_subdir}", file=sys.stderr)
' "$PACKAGE_CACHE" "$NATIVE_SUBDIR"
}
new_release_id() {
    printf '%s-%s\n' "$(date -u +%Y%m%dT%H%M%SZ)" "$$"
}

create_control_from_seed() {
    local target="$1" role="$2" specification
    local -a format_args=()
    specification="$(extract_env_spec "$role")" || return $?
    echo "使用 Seed Conda 创建 control candidate：$target"
    echo "  seed         : $SEED_CONDA"
    echo "  specification: $specification"

    if run_conda_probe "$SEED_CONDA" env create --help 2>&1 | grep -q -- '--format'; then
        format_args=(--format environment-yaml)
    fi

    env \
        "${CLEAN_ENVIRONMENT[@]}" \
        CONDA_ALWAYS_YES=true \
        CONDA_DEFAULT_THREADS=1 \
        CONDA_EXECUTE_THREADS=1 \
        CONDA_FETCH_THREADS=1 \
        CONDA_REMOTE_MAX_RETRIES=5 \
        CONDA_REPODATA_THREADS=1 \
        CONDA_VERIFY_THREADS=1 \
        PATH="$SEED_CONDA_DIR:$HOST_PATH" \
        PYTHONNOUSERSITE=1 \
        "$SEED_CONDA" env create \
            "${format_args[@]}" \
            --prefix "$target" \
            --file "$specification"
}

create_from_control() {
    local target="$1" role="$2" specification
    local -a format_args=()
    specification="$(extract_env_spec "$role")" || return $?
    echo "使用新 control Conda 创建 candidate：$target"
    echo "  control conda: $CONTROL_CONDA"
    echo "  specification: $specification"

    if run_conda_probe "$CONTROL_CONDA" env create --help 2>&1 | grep -q -- '--format'; then
        format_args=(--format environment-yaml)
    fi

    env \
        "${CLEAN_ENVIRONMENT[@]}" \
        CONDA_ALWAYS_YES=true \
        CONDA_DEFAULT_THREADS=1 \
        CONDA_EXECUTE_THREADS=1 \
        CONDA_FETCH_THREADS=1 \
        CONDA_REMOTE_MAX_RETRIES=5 \
        CONDA_REPODATA_THREADS=1 \
        CONDA_VERIFY_THREADS=1 \
        PATH="$CONTROL_CONDA_DIR:$HOST_PATH" \
        PYTHONNOUSERSITE=1 \
        "$CONTROL_CONDA" env create \
            --solver libmamba \
            "${format_args[@]}" \
            --prefix "$target" \
            --file "$specification"
}
validate_control_conda() {
    local prefix="$1" conda_bin="$1/bin/conda" version_output major
    [[ -x "$conda_bin" ]] || return 1
    version_output="$(run_in_sanitized_environment "$prefix/bin" "$conda_bin" --version 2>&1)" || return 1
    major="$(printf '%s\n' "$version_output" | sed -nE 's/^conda[[:space:]]+([0-9]+)(\..*)?$/\1/p' | head -n1)"
    [[ "$major" =~ ^[0-9]+$ ]] && (( major >= 26 )) || return 1
}

smoke_control() {
    local prefix="$1" version major
    prefix_has_tools "$prefix" "${CONTROL_REQUIRED[@]}" || return 1
    prefix_matches_native_platform "$prefix" || return 1
    validate_control_conda "$prefix" || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/python" - <<'PY' >/dev/null
from importlib.metadata import version
import jsonschema
import yaml
version("snakemake")
version("snakemake-executor-plugin-slurm")
PY
    version="$(run_in_sanitized_environment "$prefix/bin" "$prefix/bin/snakemake" --version | head -n1)" || return 1
    major="$(printf '%s' "$version" | sed -E 's/^([0-9]+).*/\1/')"
    [[ "$major" =~ ^[0-9]+$ ]] && (( major >= 9 )) || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/cargo" --version >/dev/null 2>&1 || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/rustc" --version >/dev/null 2>&1 || return 1
}

smoke_biscuit() {
    local prefix="$1"
    prefix_has_tools "$prefix" "${BISCUIT_REQUIRED[@]}" || return 1
    prefix_matches_native_platform "$prefix" || return 1
    prefix_python_imports "$prefix" yaml || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/biscuit" version >/dev/null 2>&1 || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/dupsifter" --version >/dev/null 2>&1 || return 1
}

smoke_bismark() {
    local prefix="$1"
    prefix_has_tools "$prefix" "${BISMARK_REQUIRED[@]}" || return 1
    prefix_matches_native_platform "$prefix" || return 1
    prefix_python_imports "$prefix" yaml pysam umi_tools || return 1
    [[ "$(run_in_sanitized_environment "$prefix/bin" "$prefix/bin/umi_tools" --version)" == *"1.1.6"* ]] || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/bismark" --version >/dev/null 2>&1 || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/bismark" prepare --help >/dev/null 2>&1 || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/bismark_genome_preparation" --help >/dev/null 2>&1 || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/bowtie2" --version >/dev/null 2>&1 || return 1
}

# 验证 TAPS 工具链与固定 Rastair 版本。
smoke_rastair() {
    local prefix="$1"
    prefix_has_tools "$prefix" "${RASTAIR_REQUIRED[@]}" || return 1
    prefix_matches_native_platform "$prefix" || return 1
    prefix_python_imports "$prefix" yaml pysam || return 1
    [[ "$(run_in_sanitized_environment "$prefix/bin" "$prefix/bin/rastair" --version)" == "rastair 2.2.0" ]] || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/rastair" call --help >/dev/null || return 1
}

smoke_notebook() {
    local prefix="$1"
    prefix_has_tools "$prefix" "${NOTEBOOK_REQUIRED[@]}" || return 1
    prefix_matches_native_platform "$prefix" || return 1
    prefix_python_imports "$prefix" \
        anndata IPython ipykernel jupyterlab matplotlib numba numpy packaging \
        pandas polars sklearn seaborn snapatac2 tqdm yaml || return 1
    run_in_sanitized_environment "$prefix/bin" "$prefix/bin/jupyter" --version >/dev/null || return 1
}

retry_transient_create() {
    local label="$1" target="$2"; shift 2
    local attempt=1 max_attempts=3
    while true; do
        if "$@"; then
            return 0
        fi
        if (( attempt >= max_attempts )); then
            return 1
        fi
        echo "WARNING: $label Conda create 失败；清理未完成 prefix 后重试 ($attempt/$max_attempts)..." >&2
        rm -rf -- "$target"
        sleep $(( attempt * 3 ))
        attempt=$((attempt + 1))
    done
}
create_or_reuse_role() {
    local role="$1" target="$2" smoke_fn="$3"
    local current_prefix="$CURRENT_RELEASE/$role"

    if [[ -L "$CURRENT_RELEASE" && -x "$CURRENT_RELEASE/control/bin/python" && -d "$current_prefix" ]]; then
        if run_conda_state_cli "$CURRENT_RELEASE/control/bin/python" env-reusable \
            "$current_prefix" --envs "$ENVS_SPEC" --role "$role" >/dev/null 2>&1 \
            && "$smoke_fn" "$current_prefix" >/dev/null 2>&1; then
            local resolved_current_prefix
            resolved_current_prefix="$(cd "$current_prefix" && pwd -P)"
            ln -s "$resolved_current_prefix" "$target"
            echo "复用未变化的 $role 环境：$resolved_current_prefix"
            return 0
        fi
    fi

    echo "创建/更新 $role 环境：$target"
    retry_transient_create "$role" "$target" create_from_control "$target" "$role"
}

create_all_other_envs() {
    echo "按角色增量创建/复用 biscuit / bismark / rastair / notebook 环境。"
    create_or_reuse_role biscuit "$BISCUIT_TARGET" smoke_biscuit || return 1
    create_or_reuse_role bismark "$BISMARK_TARGET" smoke_bismark || return 1
    create_or_reuse_role rastair "$RASTAIR_TARGET" smoke_rastair || return 1
    create_or_reuse_role notebook "$NOTEBOOK_TARGET" smoke_notebook || return 1
}

smoke_all_other_envs() {
    echo "统一检查 biscuit / bismark / rastair / notebook 环境。"
    smoke_biscuit "$BISCUIT_TARGET" || { echo "ERROR: biscuit 环境未通过 smoke test：$BISCUIT_TARGET" >&2; return 1; }
    smoke_bismark "$BISMARK_TARGET" || { echo "ERROR: bismark 环境未通过 smoke test：$BISMARK_TARGET" >&2; return 1; }
    smoke_rastair "$RASTAIR_TARGET" || { echo "ERROR: rastair 环境未通过 smoke test： $RASTAIR_TARGET" >&2; return 1; }
    smoke_notebook "$NOTEBOOK_TARGET" || { echo "ERROR: notebook 环境未通过 smoke test：$NOTEBOOK_TARGET" >&2; return 1; }
}

NEW_RELEASE_PATH=""
CONTROL_TARGET=""
BISCUIT_TARGET=""
BISMARK_TARGET=""
RASTAIR_TARGET=""
NOTEBOOK_TARGET=""
CONTROL_CONDA=""
CONTROL_CONDA_DIR=""
CONDA_PYTHON_BIN=""

# 从已分配的 release 根一次性派生全部 Conda 事务路径。
set_conda_release_targets() {
    [[ -n "${NEW_RELEASE_PATH:-}" ]] || {
        echo "ERROR: Conda release 路径尚未分配。" >&2
        return 2
    }
    case "$NEW_RELEASE_PATH" in
        "$RELEASE_STORE/"*) ;;
        *)
            echo "ERROR: Conda release 路径不在受管目录内：$NEW_RELEASE_PATH" >&2
            return 2
            ;;
    esac
    CONTROL_TARGET="$NEW_RELEASE_PATH/control"
    BISCUIT_TARGET="$NEW_RELEASE_PATH/biscuit"
    BISMARK_TARGET="$NEW_RELEASE_PATH/bismark"
    RASTAIR_TARGET="$NEW_RELEASE_PATH/rastair"
    NOTEBOOK_TARGET="$NEW_RELEASE_PATH/notebook"
    CONTROL_CONDA="$CONTROL_TARGET/bin/conda"
    CONTROL_CONDA_DIR="$CONTROL_TARGET/bin"
    CONDA_PYTHON_BIN="$CONTROL_TARGET/bin/python"
}

current_points_to_new_release() {
    [[ -n "${NEW_RELEASE_PATH:-}" ]] || return 1
    [[ -L "$CURRENT_RELEASE" ]] || return 1
    [[ "$(readlink "$CURRENT_RELEASE" 2>/dev/null || true)" == "$NEW_RELEASE_PATH" ]]
}

cleanup_setup_transaction() {
    local exit_status=$?
    trap - EXIT

    if [[ -n "${NEW_RELEASE_PATH:-}" ]] && ! current_points_to_new_release; then
        case "$NEW_RELEASE_PATH" in
            "$RELEASE_STORE/"*"-$$")
                if [[ -d "$NEW_RELEASE_PATH" ]]; then
                    echo "删除未发布 release：$NEW_RELEASE_PATH" >&2
                    rm -rf -- "$NEW_RELEASE_PATH" || true
                fi
                ;;
            *)
                echo "WARNING: 拒绝删除非本事务 release：$NEW_RELEASE_PATH" >&2
                ;;
        esac
    fi

    txn_release_lock "$ENVIRONMENT_LOCK" || true
    exit "$exit_status"
}

# 经指定 Python 调用 identity 模块的环境身份 / release 状态 CLI（保持净化环境）。
run_conda_state_cli() {
    local python_bin="$1"
    shift
    env \
        -u PYTHONHOME -u PYTHONSTARTUP -u PYTHONUSERBASE -u VIRTUAL_ENV \
        PYTHONPATH="$PIPELINE_ROOT/core" \
        PYTHONNOUSERSITE=1 \
        "$python_bin" -m dna_pipeline "$@"
}

check_current_release() {
    local control="$CURRENT_RELEASE/control"
    local biscuit="$CURRENT_RELEASE/biscuit"
    local bismark="$CURRENT_RELEASE/bismark"
    local notebook="$CURRENT_RELEASE/notebook"
    local rastair="$CURRENT_RELEASE/rastair"
    local python_bin="$control/bin/python"

    [[ -L "$CURRENT_RELEASE" ]] || { echo "ERROR: conda/current 缺失或不是 symlink。" >&2; return 1; }
    [[ -x "$python_bin" ]] || { echo "ERROR: current/control Python 缺失。" >&2; return 1; }

    run_conda_state_cli "$python_bin" env-validate \
        "$CURRENT_RELEASE" "$NATIVE_SUBDIR" --envs "$ENVS_SPEC" >/dev/null || return 1

    smoke_control "$control" || { echo "ERROR: current/control 未通过 smoke test。" >&2; return 1; }
    smoke_biscuit "$biscuit" || { echo "ERROR: current/biscuit 未通过 smoke test。" >&2; return 1; }
    smoke_bismark "$bismark" || { echo "ERROR: current/bismark 未通过 smoke test。" >&2; return 1; }
    smoke_rastair "$rastair" || { echo "ERROR: current/rastair 未通过 smoke test。" >&2; return 1; }
    smoke_notebook "$notebook" || { echo "ERROR: current/notebook 未通过 smoke test。" >&2; return 1; }
    echo "READY: Conda current release 已通过 $NATIVE_SUBDIR 校验。"
}

record_new_environment() {
    local prefix="$1" role="$2"
    echo "记录 $role thin environment identity。"
    run_conda_state_cli "$CONDA_PYTHON_BIN" env-record \
        "$prefix" --envs "$ENVS_SPEC" --role "$role" >/dev/null
}


finalize_release() {
    run_conda_state_cli "$CONDA_PYTHON_BIN" env-finalize \
        "$NEW_RELEASE_PATH" "$RELEASE_ID" "$NATIVE_SUBDIR" >/dev/null
}

publish_release() {
    run_conda_state_cli "$CONDA_PYTHON_BIN" env-publish \
        "$CURRENT_RELEASE" "$NEW_RELEASE_PATH" >/dev/null
}


# 保留当前、上一代 release 及其复用角色的实际目录，避免清理后产生悬空链接。
prune_superseded_releases() {
    local candidate remaining release prefix resolved referenced scan_index=0
    local -a kept=("$NEW_RELEASE_PATH")
    remaining="$(find "$RELEASE_STORE" -mindepth 1 -maxdepth 1 -type d -name '20*' \
        ! -path "$NEW_RELEASE_PATH" | LC_ALL=C sort -r)"
    while IFS= read -r candidate; do
        [[ -n "$candidate" ]] || continue
        kept+=("$candidate")
        break
    done <<< "$remaining"
    while (( scan_index < ${#kept[@]} )); do
        release="${kept[$scan_index]}"
        scan_index=$((scan_index + 1))
        for prefix in "$release"/*; do
            [[ -L "$prefix" && -d "$prefix" ]] || continue
            resolved="$(cd "$prefix" && pwd -P)" || return $?
            case "$resolved" in
                "$RELEASE_STORE"/*/*)
                    resolved="${resolved%/*}"
                    referenced=0
                    for candidate in "${kept[@]}"; do
                        [[ "$candidate" != "$resolved" ]] || { referenced=1; break; }
                    done
                    (( referenced == 1 )) || kept+=("$resolved")
                    ;;
            esac
        done
    done
    while IFS= read -r candidate; do
        [[ -n "$candidate" ]] || continue
        referenced=0
        for release in "${kept[@]}"; do
            [[ "$candidate" != "$release" ]] || { referenced=1; break; }
        done
        (( referenced == 0 )) || continue
        echo "回收过期 Conda release：$candidate"
        rm -rf -- "$candidate" || true
    done <<< "$remaining"
}



# conda 段入口：分派的唯一调用目标；flow 首行自带 set 严格模式。
run_conda_section() {
set -euo pipefail
CHECK_ONLY=0
NEW_RELEASE_PATH=""
CONTROL_TARGET=""
BISCUIT_TARGET=""
BISMARK_TARGET=""
RASTAIR_TARGET=""
NOTEBOOK_TARGET=""
CONTROL_CONDA=""
CONTROL_CONDA_DIR=""
CONDA_PYTHON_BIN=""
while (( $# > 0 )); do
    case "$1" in
        --check)
            CHECK_ONLY=1
            ;;
        -h|--help|help)
            usage_conda
            return 0
            ;;
        *)
            echo "ERROR: core/doctor.sh conda 不接受参数：$1" >&2
            echo "Seed Conda 路径由脚本自动寻找；跨组件检查/修复请使用 core/doctor.sh。" >&2
            usage_conda >&2
            return 2
            ;;
    esac
    shift
done
if [[ "$(uname -s)" == Linux && -n "${SLURM_JOB_ID:-}" && "$CHECK_ONLY" != 1 ]]; then
    echo "ERROR: core/doctor.sh conda 构建/发布不应在 SLURM 计算作业中运行；只允许 --check。" >&2
    echo "请退出当前计算节点/SLURM job，在 HPC 登录节点重新执行：" >&2
    echo "  bash core/doctor.sh conda" >&2
    return 2
fi
unset PYTHONHOME PYTHONPATH PYTHONSTARTUP PYTHONUSERBASE VIRTUAL_ENV || true
unset LD_LIBRARY_PATH LD_PRELOAD LIBRARY_PATH PKG_CONFIG_PATH || true
unset DYLD_LIBRARY_PATH CMAKE_PREFIX_PATH || true
export PYTHONNOUSERSITE=1
require_specification "$ENVS_SPEC"
if (( CHECK_ONLY == 1 )); then
    if check_current_release; then
        return 0
    else
        rc=$?
        exit "$rc"
    fi
fi
mkdir -p \
    "$CONDA_ROOT" \
    "$PACKAGE_CACHE" \
    "$SETUP_HOME" \
    "$SETUP_ENVS" \
    "$SETUP_CACHE" \
    "$SETUP_CONFIG"
SEED_CONDA="$(resolve_seed_conda)"
SEED_CONDA_DIR="$(cd "$(dirname "$SEED_CONDA")" && pwd -P)"
export PATH="$SEED_CONDA_DIR:$HOST_PATH"
validate_seed_context
txn_acquire_lock "$ENVIRONMENT_LOCK" "setup_conda transaction"
trap cleanup_setup_transaction EXIT
if [[ -e "$CURRENT_RELEASE" && ! -L "$CURRENT_RELEASE" ]]; then
    echo "ERROR: conda/current 存在但不是 symlink，拒绝继续。" >&2
    return 1
fi
mkdir -p "$RELEASE_STORE"
RELEASE_ID="$(new_release_id)"
NEW_RELEASE_PATH="$RELEASE_STORE/$RELEASE_ID"
set_conda_release_targets
if [[ -e "$NEW_RELEASE_PATH" || -L "$NEW_RELEASE_PATH" ]]; then
    echo "ERROR: release 路径已经存在：$NEW_RELEASE_PATH" >&2
    return 1
fi
mkdir "$NEW_RELEASE_PATH"
echo "创建增量 release：$NEW_RELEASE_PATH"
if [[ -L "$CURRENT_RELEASE" && -x "$CURRENT_RELEASE/control/bin/python" ]] \
    && run_conda_state_cli "$CURRENT_RELEASE/control/bin/python" env-reusable \
        "$CURRENT_RELEASE/control" --envs "$ENVS_SPEC" --role control >/dev/null 2>&1 \
    && smoke_control "$CURRENT_RELEASE/control" >/dev/null 2>&1; then
    CURRENT_CONTROL_RESOLVED="$(cd "$CURRENT_RELEASE/control" && pwd -P)"
    ln -s "$CURRENT_CONTROL_RESOLVED" "$CONTROL_TARGET"
    echo "复用未变化的 control 环境：$CURRENT_CONTROL_RESOLVED"
else
    echo "创建/更新 control：$CONTROL_TARGET"
    retry_transient_create control "$CONTROL_TARGET" create_control_from_seed "$CONTROL_TARGET" control
    smoke_control "$CONTROL_TARGET" || {
        echo "ERROR: 新 control 环境未通过 smoke test：$CONTROL_TARGET" >&2
        return 1
    }
fi
echo "control smoke test：PASS"
echo "正式环境管理器已切换为：$CONTROL_CONDA"
create_all_other_envs || {
    echo "ERROR: 其它 Conda 环境未全部创建成功。" >&2
    return 1
}
echo "其它 Conda 环境统一创建：PASS"
smoke_all_other_envs || exit 1
echo "其它 Conda 环境统一检查：PASS"
[[ -L "$CONTROL_TARGET" ]] || record_new_environment "$CONTROL_TARGET" control
[[ -L "$BISCUIT_TARGET" ]] || record_new_environment "$BISCUIT_TARGET" biscuit
[[ -L "$BISMARK_TARGET" ]] || record_new_environment "$BISMARK_TARGET" bismark
[[ -L "$RASTAIR_TARGET" ]] || record_new_environment "$RASTAIR_TARGET" rastair
[[ -L "$NOTEBOOK_TARGET" ]] || record_new_environment "$NOTEBOOK_TARGET" notebook
finalize_release
publish_release
echo "已原子发布 Conda release：$NEW_RELEASE_PATH"
prune_superseded_releases
txn_release_lock "$ENVIRONMENT_LOCK"
trap - EXIT
echo
echo "Alopex Conda release 构建并激活成功："
echo "  platform     : $NATIVE_SUBDIR"
echo "  seed         : $SEED_CONDA"
echo "  control conda: $CONTROL_CONDA"
echo "  release      : $NEW_RELEASE_PATH"
echo "  current      : $CURRENT_RELEASE -> $NEW_RELEASE_PATH"
echo "  control      : $CONTROL_TARGET"
echo "  biscuit      : $BISCUIT_TARGET"
echo "  bismark      : $BISMARK_TARGET"
echo "  rastair      : $RASTAIR_TARGET"
echo "  notebook     : $NOTEBOOK_TARGET"
echo
echo "后续状态检查 / 自动修复请统一执行："
echo "  bash $PIPELINE_ROOT/core/doctor.sh"
}






CONTROL_ENV="$CURRENT_RELEASE/control"
BISCUIT_ENV="$CURRENT_RELEASE/biscuit"
BISMARK_ENV="$CURRENT_RELEASE/bismark"
RASTAIR_ENV="$CURRENT_RELEASE/rastair"

PYTHON_BIN="$CONTROL_ENV/bin/python"
SAMTOOLS_BIN="$BISCUIT_ENV/bin/samtools"
BISCUIT_BIN="$BISCUIT_ENV/bin/biscuit"
BWA_BIN="$RASTAIR_ENV/bin/bwa"
BISMARK_BIN="$BISMARK_ENV/bin/bismark"
BISMARK_GENOME_PREP_BIN="$BISMARK_ENV/bin/bismark_genome_preparation"
BOWTIE2_BUILD_BIN="$BISMARK_ENV/bin/bowtie2-build"

SOURCE_DIR="$PIPELINE_ROOT/resources/reference_source"
HG38_SOURCE="$SOURCE_DIR/hg38.fa"
MM10_SOURCE="$SOURCE_DIR/mm10.fa"
LAMBDA_SOURCE="$SOURCE_DIR/lambda.fa"
PUC19_SOURCE="$SOURCE_DIR/pUC19.fa"
HG38_GTF_SOURCE="$SOURCE_DIR/hg38.genes.gtf"
MM10_GTF_SOURCE="$SOURCE_DIR/mm10.genes.gtf"

HG38_REFERENCE="$PIPELINE_ROOT/resources/hg38_reference/hg38.primary.lambda.puc19.fa"
MM10_REFERENCE="$PIPELINE_ROOT/resources/mm10_reference/mm10.primary.lambda.puc19.fa"
HG38_TSS_BED="$PIPELINE_ROOT/resources/hg38_reference/tss/hg38_TSS_2000_2000_20.bed"
HG38_SINGLE_CPG="$PIPELINE_ROOT/resources/hg38_reference/cpg/hg38.single_cpg.bed.gz"
MM10_TSS_BED="$PIPELINE_ROOT/resources/mm10_reference/tss/mm10_TSS_2000_2000_20.bed"
MM10_SINGLE_CPG="$PIPELINE_ROOT/resources/mm10_reference/cpg/mm10.single_cpg.bed.gz"


REFERENCE_MIN_CPUS=16
REFERENCE_MIN_MEM_MIB=$((128 * 1024))
REFERENCE_SUBMIT_CPUS=16
REFERENCE_SUBMIT_MEM="128G"
REFERENCE_LOG_ROOT="$PIPELINE_ROOT/logs"
REFERENCE_LOCK="$PIPELINE_ROOT/resources/.locks/reference-build.lock.d"
REFERENCE_BUILD_ROOT="$PIPELINE_ROOT/resources/.build"


REFERENCE_CHILD_PIDS=""

REFERENCE_TEMP_PATHS=()

register_reference_temp() {
    local path="$1"
    case "$path" in
        "$REFERENCE_BUILD_ROOT"/*) ;;
        *) echo "ERROR: 拒绝登记 $REFERENCE_BUILD_ROOT 之外的临时路径： $path" >&2; return 2 ;;
    esac
    REFERENCE_TEMP_PATHS+=("$path")
}

cleanup_reference_temps() {
    local path
    for path in ${REFERENCE_TEMP_PATHS[@]+"${REFERENCE_TEMP_PATHS[@]}"}; do
        case "$path" in
            "$REFERENCE_BUILD_ROOT"/*) rm -rf -- "$path" ;;
        esac
    done
    REFERENCE_TEMP_PATHS=()
}

stop_reference_children() {
    local child_pid
    for child_pid in ${REFERENCE_CHILD_PIDS:-}; do
        [[ "$child_pid" =~ ^[0-9]+$ ]] || continue
        kill -TERM "$child_pid" 2>/dev/null || true
    done
    for child_pid in ${REFERENCE_CHILD_PIDS:-}; do
        [[ "$child_pid" =~ ^[0-9]+$ ]] || continue
        wait "$child_pid" 2>/dev/null || true
    done
    REFERENCE_CHILD_PIDS=""
}

# 持锁且无活动构建子进程时回收候选目录，保留提交 journal 与恢复备份。
cleanup_stale_reference_builds() {
    mkdir -p "$REFERENCE_BUILD_ROOT"
    find "$REFERENCE_BUILD_ROOT" -mindepth 1 -maxdepth 1 -type d -name 'reference.*' -exec rm -rf {} + 2>/dev/null || true
}

reference_transaction_exit() {
    stop_reference_children
    recover_downstream_transactions || true
    recover_biscuit_commit || true
    recover_bwa_commit || true
    recover_bismark_commit || true
    cleanup_reference_temps
    cleanup_stale_reference_builds
    txn_release_lock "$REFERENCE_LOCK"
}

reference_transaction_signal() {
    local exit_code="$1"
    stop_reference_children
    exit "$exit_code"
}

with_transaction_lock_reference() {
    local rc
    txn_acquire_lock "$REFERENCE_LOCK" "reference build transaction" || return $?
    cleanup_stale_reference_builds
    recover_downstream_transactions && recover_biscuit_commit && recover_bismark_commit && recover_bwa_commit || {
        txn_release_lock "$REFERENCE_LOCK"
        return 1
    }
    trap 'reference_transaction_exit' EXIT
    trap 'reference_transaction_signal 130' INT
    trap 'reference_transaction_signal 143' TERM
    trap 'reference_transaction_signal 129' HUP

    if "$@"; then rc=0; else rc=$?; fi
    stop_reference_children
    cleanup_reference_temps
    cleanup_stale_reference_builds
    trap - EXIT INT TERM HUP
    txn_release_lock "$REFERENCE_LOCK"
    return "$rc"
}

REFERENCES_READY=0
HG38_ASSEMBLED_READY=0
HG38_BWA_READY=0
HG38_BISCUIT_READY=0
HG38_BISMARK_STANDARD_READY=0
HG38_BISMARK_COMBINED_READY=0
HG38_BISMARK_READY=0
HG38_DOWNSTREAM_READY=0
HG38_BUNDLE_READY=0
MM10_ASSEMBLED_READY=0
MM10_BWA_READY=0
MM10_BISCUIT_READY=0
MM10_BISMARK_STANDARD_READY=0
MM10_BISMARK_COMBINED_READY=0
MM10_BISMARK_READY=0
MM10_DOWNSTREAM_READY=0
MM10_BUNDLE_READY=0

# BISCUIT index 后缀通过 dna_pipeline index-suffixes CLI 惰性读取，避免 Shell 双写清单。
_biscuit_index_suffixes() {
    printf '.fai\n'
    env PYTHONPATH="$PIPELINE_ROOT/core" "$PYTHON_BIN" -m dna_pipeline index-suffixes
}
# index 中断时 BISCUIT 可能留下的工作文件；不属于 READY identity，只用于清理。
_biscuit_transient_suffixes() {
    env PYTHONPATH="$PIPELINE_ROOT/core" "$PYTHON_BIN" -m dna_pipeline index-suffixes --transient
}

usage_reference() {
    cat <<'USAGE'
用法：
  bash core/doctor.sh reference [--check]

作用：
  构建、验证并记录 Alopex 标准 hg38 / mm10 reference。

参数：
  --check
      只运行轻量 preflight 并验证两套标准 bundle 是否 READY；不构建、不提交 SLURM。

仅在需要重建对应组件时要求 source：
  assembled FASTA rebuild: 对应 hg38/mm10 genome + lambda + pUC19
  downstream rebuild: 对应 species 的固定 GENCODE GTF
  现成完整 FASTA/index 可直接导入并复用，无需 source 文件。

标准输出：
  hg38: chr1-22,chrX,chrY,chrM + lambda + pUC19
  mm10: chr1-19,chrX,chrY,chrM + lambda + pUC19
  下游: protein-coding TSS ±2000 bp / 20 bp bins + primary-contig single-CpG BED.gz

本脚本不会联网下载 reference。
其它物种由用户自行准备最终 FASTA，并自行构建 FAI、BISCUIT、Bismark index 与下游 reference。

运行前检查：
  软件/配置 + 4 个 source FASTA / 2 个 GTF + 已有 reference bundle 是否全部 READY。
  source FASTA 通过 .fai 检查必需 contig 与预期长度；不会为完整性检查全量哈希序列。
  preflight 不生成详细 build plan；若 hg38/mm10 已全部 READY，脚本直接成功结束，不申请 SLURM。

Linux/HPC 行为：
  当前有效 SLURM allocation 只要 >= 16 CPU / 128 GiB 就直接构建；
  因此 32 CPU / 128 GiB、16 CPU / 256 GiB 等更大 allocation 都直接复用。
  若不在 SLURM、资源不足或无法安全确认资源，则申请 16 CPU / 128 GiB，
  新 job 重新执行本脚本并再次经过同一 preflight/resource gate。
  若已有小资源 allocation 禁止 nested sbatch，请回登录节点重跑。

并发约束：
  同一 Alopex/reference 目录不支持并发执行多个 core/doctor.sh reference。
  请保证一次只运行一个 reference setup；违反该约束属于不支持的用法。
  单次受锁事务会按 CPU/RAM 自适应并行两个物种的独立阶段。

失败清理：
  所有 rebuild 都在 resources/.build/reference.* 临时目录完成并验证后提交。
  FAIL、Ctrl-C、TERM 或 HUP 只清理未提交临时目录；运行前已存在的正式 reference/index 不作为 rollback workspace。
USAGE
}



check_prerequisites_reference() {
    local mode="${1:-readonly}"
    [[ -x "$PYTHON_BIN" ]] || {
        echo "ERROR: control Python 缺失：$PYTHON_BIN" >&2
        echo "请先运行：bash core/doctor.sh" >&2
        return 127
    }
    [[ "$mode" == readonly ]] && return 0
    [[ "$mode" == build ]] || { echo "ERROR: 无效的前置检查模式： $mode" >&2; return 2; }
    [[ -x "$BWA_BIN" ]] || { echo "ERROR: reference 构建需要 rastair 环境的 bwa： $BWA_BIN" >&2; return 127; }
    [[ -x "$SAMTOOLS_BIN" && -x "$BISCUIT_BIN" ]] || {
        echo "ERROR: reference 构建需要 biscuit 环境的 samtools/biscuit： $BISCUIT_ENV" >&2
        return 127
    }
    [[ -x "$BISMARK_BIN" && -x "$BISMARK_GENOME_PREP_BIN" && -x "$BOWTIE2_BUILD_BIN" ]] || {
        echo "ERROR: reference 构建需要 bismark/bismark_genome_preparation/bowtie2-build： $BISMARK_ENV" >&2
        return 127
    }
    mkdir -p "$SOURCE_DIR" \
        "$(dirname "$HG38_REFERENCE")" "$(dirname "$MM10_REFERENCE")" \
        "$(dirname "$HG38_TSS_BED")" "$(dirname "$MM10_TSS_BED")" \
        "$(dirname "$HG38_SINGLE_CPG")" "$(dirname "$MM10_SINGLE_CPG")"
}


fai_contract_check() {
    local mode="$1" label="$2" fai="$3" fasta="$4"
    env PYTHONPATH="$PIPELINE_ROOT/core" "$PYTHON_BIN" -m dna_pipeline fai-check \
        --mode "$mode" --label "$label" --fai "$fai" --fasta "$fasta" >/dev/null 2>&1
}

ensure_source_fai() {
    local label="$1" fasta="$2" fai="${2}.fai" rebuilt=0

    if [[ ! -s "$fai" || "$fasta" -nt "$fai" ]]; then
        rm -f -- "$fai"
        echo "为 source FASTA 建索引以检查完整性： $fasta"
        "$SAMTOOLS_BIN" faidx "$fasta"
        rebuilt=1
    fi

    if fai_contract_check source "$label" "$fai" "$fasta"; then
        return 0
    fi

    if (( rebuilt == 0 )); then
        echo "缓存的 source FAI 校验失败，重建一次： $fai"
        rm -f -- "$fai"
        "$SAMTOOLS_BIN" faidx "$fasta"
        fai_contract_check source "$label" "$fai" "$fasta" && return 0
    fi

    echo "ERROR: source FASTA 未通过轻量完整性检查：$fasta" >&2
    echo "       必需 contig 缺失或长度不符合标准；文件可能下载截断或不是预期 reference。" >&2
    return 2
}

check_sources_for_build_plan() {
    local missing=0

    if (( HG38_ASSEMBLED_READY == 0 )); then
        for path in "$HG38_SOURCE" "$LAMBDA_SOURCE" "$PUC19_SOURCE"; do
            [[ -s "$path" ]] || { echo "缺少或为空：$path" >&2; missing=1; }
        done
        (( missing == 0 )) || return 2
        ensure_source_fai hg38 "$HG38_SOURCE" || return $?
        ensure_source_fai lambda "$LAMBDA_SOURCE" || return $?
        ensure_source_fai pUC19 "$PUC19_SOURCE" || return $?
    fi
    if (( MM10_ASSEMBLED_READY == 0 )); then
        for path in "$MM10_SOURCE" "$LAMBDA_SOURCE" "$PUC19_SOURCE"; do
            [[ -s "$path" ]] || { echo "缺少或为空：$path" >&2; missing=1; }
        done
        (( missing == 0 )) || return 2
        ensure_source_fai mm10 "$MM10_SOURCE" || return $?
        ensure_source_fai lambda "$LAMBDA_SOURCE" || return $?
        ensure_source_fai pUC19 "$PUC19_SOURCE" || return $?
    fi
    if (( HG38_DOWNSTREAM_READY == 0 )); then
        [[ -s "$HG38_GTF_SOURCE" ]] || { echo "缺少或为空：$HG38_GTF_SOURCE" >&2; return 2; }
    fi
    if (( MM10_DOWNSTREAM_READY == 0 )); then
        [[ -s "$MM10_GTF_SOURCE" ]] || { echo "缺少或为空：$MM10_GTF_SOURCE" >&2; return 2; }
    fi
    return 0
}


assembled_reference_ready() {
    local species="$1" fasta="$2"
    [[ -s "$fasta" ]] || return 1

    [[ -s "$fasta.fai" ]] || return 1
    primary_contigs_ready "$species" "$fasta"
}

primary_contigs_ready() {
    local species="$1" fasta="$2"
    [[ -s "$fasta.fai" ]] || return 1
    fai_contract_check assembled "$species" "$fasta.fai" "$fasta"
}

clear_reference_derivatives() {
    local fasta="$1" suffix
    while IFS= read -r suffix; do rm -f -- "${fasta}${suffix}"; done < <(_biscuit_index_suffixes)
    while IFS= read -r suffix; do rm -f -- "${fasta}${suffix}"; done < <(_biscuit_transient_suffixes)
    rm -rf -- "$fasta.bismark" "$fasta.bwa"
}

assemble_reference() {
    local species="$1" genome_source="$2" output="$3" build_dir tmp_fasta
    local -a wanted

    case "$species" in
        hg38) wanted=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM) ;;
        mm10) wanted=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY chrM) ;;
        *) echo "ERROR: 不支持的标准物种： $species" >&2; return 2 ;;
    esac

    mkdir -p "$REFERENCE_BUILD_ROOT" "$(dirname "$output")"
    build_dir="$(mktemp -d "$REFERENCE_BUILD_ROOT/reference.${species}.fasta.XXXXXX")" || return $?
    register_reference_temp "$build_dir" || return $?
    tmp_fasta="$build_dir/${output##*/}"
    echo "在临时目录组装 $species FASTA： $build_dir"

    "$SAMTOOLS_BIN" faidx "$genome_source" "${wanted[@]}" > "$tmp_fasta" || return $?
    env PYTHONPATH="$PIPELINE_ROOT/core" "$PYTHON_BIN" -m dna_pipeline append-spikein \
        --lambda-fasta "$LAMBDA_SOURCE" --puc19-fasta "$PUC19_SOURCE" --target-fasta "$tmp_fasta" || return $?
    "$SAMTOOLS_BIN" faidx "$tmp_fasta" || return $?
    fai_contract_check assembled "$species" "$tmp_fasta.fai" "$tmp_fasta" || {
        echo "ERROR: 临时组装的 $species FASTA 未通过 contig、顺序或长度校验。" >&2
        return 1
    }

    echo "临时组装的 $species FASTA 已通过校验，开始提交。"
    clear_reference_derivatives "$output"
    downstream_reference_paths "$species"
    rm -f -- "$DOWNSTREAM_TSS" "$DOWNSTREAM_SINGLE_CPG"
    mv -f -- "$tmp_fasta" "$output" || return $?
    mv -f -- "$tmp_fasta.fai" "$output.fai" || return $?
    rm -rf -- "$build_dir"
}


biscuit_index_ready() {
    local fasta="$1" suffix journal
    [[ -s "$fasta" ]] || return 1
    journal="$(biscuit_commit_journal "$fasta")" || return $?
    [[ ! -d "$journal" ]] || return 1
    while IFS= read -r suffix; do
        [[ -s "${fasta}${suffix}" ]] || return 1
    done < <(_biscuit_index_suffixes)
    return 0
}

# 每个物种独占 sidecar journal，允许 hg38/mm10 并行提交与独立恢复。
biscuit_commit_journal() {
    local species
    case "$1" in
        "$HG38_REFERENCE") species=hg38 ;;
        "$MM10_REFERENCE") species=mm10 ;;
        *) echo "ERROR: BISCUIT reference 不属于受管物种：$1" >&2; return 2 ;;
    esac
    printf '%s/%s.biscuit-sidecar.COMMITTING\n' "$REFERENCE_BUILD_ROOT" "$species"
}

# 恢复中断的 sidecar 提交：BACKUPS_READY 之前 public 文件从未被触碰，journal 只是 scratch。
recover_biscuit_commit() {
    local journal fasta target backup suffix
    (( $# > 0 )) || set -- "$HG38_REFERENCE" "$MM10_REFERENCE"
    for fasta in "$@"; do
        journal="$(biscuit_commit_journal "$fasta")" || return $?
        [[ -d "$journal" ]] || continue
        if [[ ! -f "$journal/BACKUPS_READY" ]]; then
            rm -rf -- "$journal" || return $?
            continue
        fi
        target="$(cat "$journal/TARGET" 2>/dev/null || true)"
        if [[ "$target" != "$fasta" ]]; then
            echo "ERROR: BISCUIT sidecar journal TARGET 与物种不符，拒绝恢复：$journal" >&2
            return 1
        fi
        echo "恢复中断的 BISCUIT sidecar 提交： $fasta"
        backup="$journal/backup"
        while IFS= read -r suffix; do
            rm -f -- "${fasta}${suffix}" || return $?
            [[ ! -e "$backup/sidecar${suffix}" ]] || cp -p -- "$backup/sidecar${suffix}" "${fasta}${suffix}" || return $?
        done < <(_biscuit_index_suffixes)
        rm -rf -- "$journal" || return $?
    done
}

commit_biscuit_sidecars() {
    local fasta="$1" tmp_fasta="$2" suffix
    while IFS= read -r suffix; do rm -f -- "${fasta}${suffix}"; done < <(_biscuit_transient_suffixes)
    while IFS= read -r suffix; do
        mv -f -- "${tmp_fasta}${suffix}" "${fasta}${suffix}" || return 1
    done < <(_biscuit_index_suffixes)
    while IFS= read -r suffix; do
        [[ -s "${fasta}${suffix}" ]] || return 1
    done < <(_biscuit_index_suffixes)
}

# 检查普通 BWA 索引的完整五文件家族。
bwa_index_ready() {
    local fasta="$1" root="${2:-$1.bwa}" path paths
    [[ "$root" != "$fasta.bwa" || ! -e "$REFERENCE_BUILD_ROOT/bwa-index.COMMITTING" ]] || return 1
    paths="$(env PYTHONPATH="$PIPELINE_ROOT/core" "$PYTHON_BIN" -m dna_pipeline index-paths --backend rastair --reference "$fasta")" || return $?
    [[ -n "$paths" ]] || return 1
    while IFS= read -r path; do
        [[ -s "$root/${path##*/}" ]] || return 1
    done <<< "$paths"
}

# 恢复未完成的普通 BWA 目录替换；完整候选提交后才删除旧索引。
recover_bwa_commit() {
    local journal="$REFERENCE_BUILD_ROOT/bwa-index.COMMITTING" target
    [[ -d "$journal" ]] || return 0
    if [[ -f "$journal/PREPARED" && ! -f "$journal/FINISHED" ]]; then
        target="$(cat "$journal/TARGET")" || return $?
        case "$target" in "$HG38_REFERENCE.bwa"|"$MM10_REFERENCE.bwa") ;; *) return 1 ;; esac
        if [[ -d "$journal/backup" || -L "$journal/backup" ]]; then
            rm -rf -- "$target" || return $?
            mv -- "$journal/backup" "$target" || return $?
        elif [[ -f "$journal/ABSENT" ]]; then
            rm -rf -- "$target" || return $?
        else
            [[ -d "$target" ]] || return 1
        fi
    fi
    rm -rf -- "$journal"
}

# 在独立目录构建普通 BWA 索引，并通过可恢复事务替换完整家族。
build_bwa_index() {
    local fasta="$1" build_dir journal="$REFERENCE_BUILD_ROOT/bwa-index.COMMITTING"
    mkdir -p "$REFERENCE_BUILD_ROOT" || return $?
    build_dir="$(mktemp -d "$REFERENCE_BUILD_ROOT/reference.bwa.XXXXXX")" || return $?
    register_reference_temp "$build_dir" || return $?
    mkdir "$build_dir/new" || return $?
    "$BWA_BIN" index -p "$build_dir/new/genome" "$fasta" || return $?
    bwa_index_ready "$fasta" "$build_dir/new" || return 1
    [[ ! -e "$journal" ]] || return 1
    mkdir "$journal" || return $?
    printf '%s\n' "$fasta.bwa" > "$journal/TARGET" || return $?
    if [[ ! -e "$fasta.bwa" && ! -L "$fasta.bwa" ]]; then : > "$journal/ABSENT"; fi
    : > "$journal/PREPARED" || return $?
    if [[ ! -f "$journal/ABSENT" ]]; then mv -- "$fasta.bwa" "$journal/backup" || return $?; fi
    mv -- "$build_dir/new" "$fasta.bwa" || return $?
    : > "$journal/FINISHED" || return $?
    rm -rf -- "$journal" "$build_dir"
}

build_biscuit_index() {
    local fasta="$1" species build_dir tmp_fasta journal backup suffix
    journal="$(biscuit_commit_journal "$fasta")" || return $?
    case "$fasta" in "$HG38_REFERENCE") species=hg38 ;; "$MM10_REFERENCE") species=mm10 ;; esac
    mkdir -p "$REFERENCE_BUILD_ROOT"
    recover_biscuit_commit "$fasta" || return $?
    build_dir="$(mktemp -d "$REFERENCE_BUILD_ROOT/reference.${species}.biscuit.XXXXXX")" || return $?
    register_reference_temp "$build_dir" || return $?
    tmp_fasta="$build_dir/${fasta##*/}"
    echo "在临时目录构建 BISCUIT index： $build_dir"
    cp "$fasta" "$tmp_fasta" || return $?
    "$SAMTOOLS_BIN" faidx "$tmp_fasta" || return $?
    "$BISCUIT_BIN" index "$tmp_fasta" || return $?
    while IFS= read -r suffix; do
        [[ -s "${tmp_fasta}${suffix}" ]] || {
            echo "ERROR: 临时 BISCUIT index 不完整，保留已有公开 index。" >&2
            return 1
        }
    done < <(_biscuit_index_suffixes)

    echo "临时 BISCUIT index 已通过校验，开始提交 sidecars。"
    [[ ! -e "$journal" ]] || {
        echo "ERROR: BISCUIT sidecar 提交 journal 尚未恢复： $journal" >&2
        return 1
    }
    mv -- "$build_dir" "$journal" || return $?
    tmp_fasta="$journal/${fasta##*/}"
    backup="$journal/backup"
    mkdir -p "$backup" || return $?
    printf '%s\n' "$fasta" > "$journal/TARGET"
    while IFS= read -r suffix; do
        [[ ! -e "${fasta}${suffix}" ]] || cp -p -- "${fasta}${suffix}" "$backup/sidecar${suffix}" || return $?
    done < <(_biscuit_index_suffixes)
    : > "$journal/BACKUPS_READY"

    if ! commit_biscuit_sidecars "$fasta" "$tmp_fasta"; then
        echo "ERROR: BISCUIT sidecar 提交失败，正在恢复上一完整文件族。" >&2
        recover_biscuit_commit "$fasta" || return $?
        return 1
    fi
    rm -rf -- "$journal"
}


# 复用统一枚举检查当前模式的完整 Bowtie2 家族；构建期校验候选目录。
bismark_component_ready() {
    local fasta="$1" root="$2" mode="$3" fasta_size genome_size paths path
    local -a args=(--reference "$fasta" --backend bismark --bismark-root "$root")
    [[ "$root" != "$fasta.bismark" || ! -d "$(bismark_commit_journal)" ]] || return 1
    [[ -s "$root/genome.fa" ]] || return 1
    fasta_size="$(stat -c %s "$fasta" 2>/dev/null || stat -f %z "$fasta")" || return 1
    genome_size="$(stat -c %s "$root/genome.fa" 2>/dev/null || stat -f %z "$root/genome.fa")" || return 1
    [[ "$fasta_size" == "$genome_size" ]] || return 1
    [[ "$mode" != local ]] || args+=(--bismark-local-alignment)
    paths="$(PYTHONPATH="$PIPELINE_ROOT/core" "$PYTHON_BIN" -m dna_pipeline index-paths "${args[@]}")" || return $?
    [[ -n "$paths" ]] || return 1
    while IFS= read -r path; do [[ -s "$path" ]] || return 1; done <<< "$paths"
}

bismark_standard_index_ready() {
    bismark_component_ready "$1" "${2:-$1.bismark}" local
}

bismark_combined_index_ready() {
    bismark_component_ready "$1" "${2:-$1.bismark}" combined
}

bismark_index_ready() {
    local fasta="$1"
    bismark_standard_index_ready "$fasta" && bismark_combined_index_ready "$fasta"
}

# Bismark 组件只移动当前模式的文件，保留另一模式；journal 不参与临时目录回收。
bismark_commit_journal() {
    printf '%s/bismark-index.COMMITTING\n' "$REFERENCE_BUILD_ROOT"
}

# 已搬回的备份以原位文件为恢复边界；恢复再次中断时不会删除已经复原的组件。
recover_bismark_commit() {
    local journal root relative existed backup target
    journal="$(bismark_commit_journal)"
    [[ -d "$journal" ]] || return 0
    if [[ ! -f "$journal/PREPARED" || -f "$journal/FINISHED" ]]; then
        rm -rf -- "$journal"
        return $?
    fi
    root="$(cat "$journal/TARGET")" || return $?
    [[ -n "$root" && -s "$journal/ENTRIES" ]] || return 1
    while IFS=$'\t' read -r relative existed; do
        case "$relative" in genome.fa|Bisulfite_Genome/CT_conversion|Bisulfite_Genome/GA_conversion|Bisulfite_Genome/Combined) ;; *) return 1 ;; esac
        target="$root/$relative"
        backup="$journal/backup/$relative"
        if [[ "$existed" == 1 ]]; then
            if [[ -e "$backup" || -L "$backup" ]]; then
                rm -rf -- "$target" || return $?
                mkdir -p "$(dirname "$target")" || return $?
                mv -- "$backup" "$target" || return $?
            else
                [[ -e "$target" || -L "$target" ]] || return 1
            fi
        elif [[ "$existed" == 0 ]]; then
            rm -rf -- "$target" || return $?
        else
            return 1
        fi
    done < "$journal/ENTRIES"
    : > "$journal/FINISHED" || return $?
    rm -rf -- "$journal"
}

commit_bismark_component() {
    local fasta="$1" candidate="$2" mode="$3" root="$1.bismark" journal relative existed
    local -a components=(genome.fa)
    case "$mode" in
        local) components+=(Bisulfite_Genome/CT_conversion Bisulfite_Genome/GA_conversion) ;;
        combined) components+=(Bisulfite_Genome/Combined) ;;
        *) return 2 ;;
    esac
    journal="$(bismark_commit_journal)"
    [[ ! -e "$journal" ]] || return 1
    mkdir -p "$journal/backup/Bisulfite_Genome" "$root/Bisulfite_Genome" || return $?
    printf '%s\n' "$root" > "$journal/TARGET" || return $?
    for relative in "${components[@]}"; do
        [[ -e "$candidate/$relative" ]] || return 1
        existed=0
        [[ ! -e "$root/$relative" && ! -L "$root/$relative" ]] || existed=1
        printf '%s\t%s\n' "$relative" "$existed" >> "$journal/ENTRIES" || return $?
    done
    : > "$journal/PREPARED" || return $?
    for relative in "${components[@]}"; do
        if { [[ ! -e "$root/$relative" && ! -L "$root/$relative" ]] || mv -- "$root/$relative" "$journal/backup/$relative"; } &&
           mv -- "$candidate/$relative" "$root/$relative"; then
            continue
        fi
        recover_bismark_commit || return $?
        return 1
    done
    : > "$journal/FINISHED" || return $?
    rm -rf -- "$journal"
}

build_bismark_standard_index() {
    local fasta="$1" species build_dir tmp_root
    case "$fasta" in "$HG38_REFERENCE") species=hg38 ;; "$MM10_REFERENCE") species=mm10 ;; *) species=custom ;; esac
    mkdir -p "$REFERENCE_BUILD_ROOT"
    build_dir="$(mktemp -d "$REFERENCE_BUILD_ROOT/reference.${species}.bismark-standard.XXXXXX")" || return $?
    register_reference_temp "$build_dir" || return $?
    tmp_root="$build_dir/genome"
    mkdir -p "$tmp_root" || return $?
    cp "$fasta" "$tmp_root/genome.fa" || return $?
    echo "在临时目录构建标准 Bismark CT/GA Bowtie2 indexes： $tmp_root"
    echo "Bismark index 线程：CT/GA 各最多 $BISMARK_INDEX_THREADS_PER_STRAND 个 bowtie2-build 线程，总计最多 $((BISMARK_INDEX_THREADS_PER_STRAND * 2)) 个。"
    if (( BISMARK_INDEX_THREADS_PER_STRAND >= 2 )); then
        "$BISMARK_GENOME_PREP_BIN" --bowtie2 --parallel "$BISMARK_INDEX_THREADS_PER_STRAND" --path_to_aligner "$(dirname "$BOWTIE2_BUILD_BIN")" "$tmp_root" || return $?
    else
        "$BISMARK_GENOME_PREP_BIN" --bowtie2 --path_to_aligner "$(dirname "$BOWTIE2_BUILD_BIN")" "$tmp_root" || return $?
    fi
    bismark_standard_index_ready "$fasta" "$tmp_root" || {
        echo "ERROR: 临时 Bismark CT/GA indexes 不完整，保留已有公开 index。" >&2
        return 1
    }

    echo "临时标准 Bismark indexes 已通过校验，开始提交 CT/GA 组件。"
    commit_bismark_component "$fasta" "$tmp_root" local || return $?
    rm -rf -- "$build_dir"
}

build_bismark_combined_index() {
    local fasta="$1" species build_dir tmp_root
    case "$fasta" in "$HG38_REFERENCE") species=hg38 ;; "$MM10_REFERENCE") species=mm10 ;; *) species=custom ;; esac
    mkdir -p "$REFERENCE_BUILD_ROOT"
    build_dir="$(mktemp -d "$REFERENCE_BUILD_ROOT/reference.${species}.bismark-combined.XXXXXX")" || return $?
    register_reference_temp "$build_dir" || return $?
    tmp_root="$build_dir/genome"
    mkdir -p "$tmp_root" || return $?
    cp "$fasta" "$tmp_root/genome.fa" || return $?
    echo "为 Combined 发布构建临时 Bismark index bundle： $tmp_root"
    echo "Bismark Combined index 最多使用 $BISMARK_COMBINED_INDEX_THREADS 个 bowtie2-build 线程。"
    if (( BISMARK_COMBINED_INDEX_THREADS >= 2 )); then
        "$BISMARK_BIN" prepare --combined_genome --parallel "$BISMARK_COMBINED_INDEX_THREADS" --path_to_aligner "$(dirname "$BOWTIE2_BUILD_BIN")" "$tmp_root" || return $?
    else
        "$BISMARK_BIN" prepare --combined_genome --path_to_aligner "$(dirname "$BOWTIE2_BUILD_BIN")" "$tmp_root" || return $?
    fi
    bismark_combined_index_ready "$fasta" "$tmp_root" || {
        echo "ERROR: 临时 Bismark Combined index 不完整，保留已有公开 CT/GA indexes。" >&2
        return 1
    }

    commit_bismark_component "$fasta" "$tmp_root" combined || return $?
    rm -rf -- "$build_dir"
    bismark_combined_index_ready "$fasta" || {
        echo "ERROR: 已提交的 Bismark Combined index 未通过校验。" >&2
        return 1
    }
    echo "Bismark Combined index 已校验并提交；已有 CT/GA indexes 已复用。"
}


downstream_reference_paths() {
    local species="$1"
    case "$species" in
        hg38)
            DOWNSTREAM_GTF="$HG38_GTF_SOURCE"
            DOWNSTREAM_TSS="$HG38_TSS_BED"
            DOWNSTREAM_SINGLE_CPG="$HG38_SINGLE_CPG"
            ;;
        mm10)
            DOWNSTREAM_GTF="$MM10_GTF_SOURCE"
            DOWNSTREAM_TSS="$MM10_TSS_BED"
            DOWNSTREAM_SINGLE_CPG="$MM10_SINGLE_CPG"
            ;;
        *)
            echo "ERROR: 不支持的受管 reference 物种： $species" >&2
            return 2
            ;;
    esac
}

downstream_references_ready() {
    local species="$1" fasta="$2"
    downstream_reference_paths "$species"

    [[ ! -e "$(downstream_commit_journal "$species")" ]] || return 1
    downstream_files_ready "$DOWNSTREAM_TSS" "$DOWNSTREAM_SINGLE_CPG"
}

downstream_files_ready() {
    local tss="$1" single_cpg="$2"
    [[ -s "$tss" ]] || return 1
    [[ -s "$single_cpg" ]] || return 1
    gzip -t "$single_cpg" >/dev/null 2>&1 || return 1
}

downstream_commit_journal() {
    local species="$1"
    printf '%s/%s.downstream.COMMITTING\n' "$REFERENCE_BUILD_ROOT" "$species"
}

recover_downstream_commit() {
    local species="$1" journal backup
    journal="$(downstream_commit_journal "$species")"
    [[ -d "$journal" ]] || return 0
    downstream_reference_paths "$species" || return $?

    if [[ ! -f "$journal/BACKUPS_READY" ]]; then
        rm -rf -- "$journal"
        return 0
    fi

    echo "恢复中断的 $species downstream reference 提交。"
    rm -f -- "$DOWNSTREAM_TSS" "$DOWNSTREAM_SINGLE_CPG" || return $?
    mkdir -p "$(dirname "$DOWNSTREAM_TSS")" "$(dirname "$DOWNSTREAM_SINGLE_CPG")" || return $?
    backup="$journal/backup"
    [[ ! -e "$backup/tss.bed" ]] || cp -p -- "$backup/tss.bed" "$DOWNSTREAM_TSS" || return $?
    [[ ! -e "$backup/single_cpg.bed.gz" ]] || cp -p -- "$backup/single_cpg.bed.gz" "$DOWNSTREAM_SINGLE_CPG" || return $?
    rm -rf -- "$journal"
}

recover_downstream_transactions() {
    recover_downstream_commit hg38 || return $?
    recover_downstream_commit mm10 || return $?
}

build_downstream_references() {
    local species="$1" fasta="$2" build_dir journal new backup
    downstream_reference_paths "$species"
    mkdir -p "$REFERENCE_BUILD_ROOT"
    recover_downstream_commit "$species" || return $?
    build_dir="$(mktemp -d "$REFERENCE_BUILD_ROOT/reference.${species}.downstream.XXXXXX")" || return $?
    register_reference_temp "$build_dir" || return $?
    new="$build_dir/new"
    mkdir -p "$new" || return $?

    echo "在临时目录构建 $species downstream TSS + single-CpG references： $build_dir"
    env PYTHONPATH="$PIPELINE_ROOT/core" "$PYTHON_BIN" -m dna_pipeline build-downstream \
        --species "$species" \
        --fasta "$fasta" \
        --annotation-gtf "$DOWNSTREAM_GTF" \
        --tss-bed "$new/tss.bed" \
        --single-cpg-bed "$new/single_cpg.bed.gz" || return $?
    downstream_files_ready "$new/tss.bed" "$new/single_cpg.bed.gz" || {
        echo "ERROR: 临时 $species downstream references 不完整，保留已有公开文件。" >&2
        return 1
    }

    journal="$(downstream_commit_journal "$species")"
    [[ ! -e "$journal" ]] || {
        echo "ERROR: downstream 提交 journal 尚未恢复： $journal" >&2
        return 1
    }
    mv -- "$build_dir" "$journal" || return $?
    new="$journal/new"
    backup="$journal/backup"
    mkdir -p "$backup" || return $?

    [[ ! -e "$DOWNSTREAM_TSS" ]] || cp -p -- "$DOWNSTREAM_TSS" "$backup/tss.bed" || return $?
    [[ ! -e "$DOWNSTREAM_SINGLE_CPG" ]] || cp -p -- "$DOWNSTREAM_SINGLE_CPG" "$backup/single_cpg.bed.gz" || return $?
    : > "$journal/BACKUPS_READY"

    mkdir -p "$(dirname "$DOWNSTREAM_TSS")" "$(dirname "$DOWNSTREAM_SINGLE_CPG")" || return $?
    if ! mv -f -- "$new/tss.bed" "$DOWNSTREAM_TSS" ||
       ! mv -f -- "$new/single_cpg.bed.gz" "$DOWNSTREAM_SINGLE_CPG" ||
       ! downstream_files_ready "$DOWNSTREAM_TSS" "$DOWNSTREAM_SINGLE_CPG"; then
        echo "ERROR: $species downstream 提交失败，正在恢复上一完整代。" >&2
        recover_downstream_commit "$species" || return $?
        return 1
    fi

    rm -rf -- "$journal"
    echo "$species downstream TSS + single-CpG 已校验并提交。"
}


validate_reference() {
    local species="$1" fasta="$2"

    primary_contigs_ready "$species" "$fasta" || {
        echo "ERROR: $species assembled FASTA 的主染色体/spike-in 结构不正确：$fasta" >&2
        return 1
    }
    biscuit_index_ready "$fasta" || {
        echo "ERROR: $species FASTA/FAI/BISCUIT index 不完整：$fasta" >&2
        return 1
    }
    bwa_index_ready "$fasta" || { echo "ERROR: BWA index 不完整： $fasta.bwa" >&2; return 1; }
    bismark_index_ready "$fasta" || {
        echo "ERROR: $species Bismark index 不完整或 genome.fa 与 reference 不一致：$fasta.bismark" >&2
        return 1
    }
    downstream_references_ready "$species" "$fasta" || {
        echo "ERROR: $species 下游 TSS/single-CpG reference 不完整：$(dirname "$fasta")" >&2
        return 1
    }
}

reference_bundle_ready() {
    local species="$1" fasta="$2"

    assembled_reference_ready "$species" "$fasta" &&
        bwa_index_ready "$fasta" &&
        biscuit_index_ready "$fasta" &&
        bismark_index_ready "$fasta" &&
        downstream_references_ready "$species" "$fasta"
}

evaluate_reference_state() {
    local species="$1" fasta="$2" prefix
    local bwa_ready=0 assembled_ready=0 biscuit_ready=0 bismark_standard_ready=0 bismark_combined_ready=0 bismark_ready=0 downstream_ready=0 bundle_ready=0

    case "$species" in
        hg38) prefix="HG38" ;;
        mm10) prefix="MM10" ;;
        *)
            echo "ERROR: 不支持的受管 reference 物种： $species" >&2
            return 2
            ;;
    esac

    assembled_reference_ready "$species" "$fasta" && assembled_ready=1 || true
    bwa_index_ready "$fasta" && bwa_ready=1 || true
    biscuit_index_ready "$fasta" && biscuit_ready=1 || true
    bismark_standard_index_ready "$fasta" && bismark_standard_ready=1 || true
    bismark_combined_index_ready "$fasta" && bismark_combined_ready=1 || true
    if (( bismark_standard_ready == 1 && bismark_combined_ready == 1 )); then bismark_ready=1; fi
    downstream_references_ready "$species" "$fasta" && downstream_ready=1 || true

    if (( assembled_ready == 1 && bwa_ready == 1 && biscuit_ready == 1 && bismark_ready == 1 && downstream_ready == 1 )); then
        bundle_ready=1
    fi

    printf -v "${prefix}_BWA_READY" '%d' "$bwa_ready"
    printf -v "${prefix}_ASSEMBLED_READY" '%d' "$assembled_ready"
    printf -v "${prefix}_BISCUIT_READY" '%d' "$biscuit_ready"
    printf -v "${prefix}_BISMARK_STANDARD_READY" '%d' "$bismark_standard_ready"
    printf -v "${prefix}_BISMARK_COMBINED_READY" '%d' "$bismark_combined_ready"
    printf -v "${prefix}_BISMARK_READY" '%d' "$bismark_ready"
    printf -v "${prefix}_DOWNSTREAM_READY" '%d' "$downstream_ready"
    printf -v "${prefix}_BUNDLE_READY" '%d' "$bundle_ready"
}

prepare_reference_build_plan() {
    evaluate_reference_state hg38 "$HG38_REFERENCE"
    evaluate_reference_state mm10 "$MM10_REFERENCE"

    if (( HG38_BUNDLE_READY == 1 && MM10_BUNDLE_READY == 1 )); then
        echo "构建计划：hg38 / mm10 reference bundles 已 READY。"
    else
        echo "构建计划：仅修复未通过校验的实际组件。"
    fi
}
detect_reference_build_resources() {
    local platform cpus="" mem_mib="" mem_bytes="" detected
    platform="$(uname -s)"

    if [[ "${SLURM_CPUS_ON_NODE:-}" =~ ^[1-9][0-9]*$ ]]; then
        cpus="$SLURM_CPUS_ON_NODE"
    elif [[ "${SLURM_CPUS_PER_TASK:-}" =~ ^[1-9][0-9]*$ ]]; then
        cpus="$SLURM_CPUS_PER_TASK"
    elif [[ "$platform" == Darwin ]]; then
        cpus="$(/usr/sbin/sysctl -n hw.ncpu 2>/dev/null || true)"
    else
        cpus="$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)"
    fi
    [[ "$cpus" =~ ^[1-9][0-9]*$ ]] || cpus=1
    (( cpus <= 16 )) || cpus=16

    if [[ "${SLURM_MEM_PER_NODE:-}" =~ ^[1-9][0-9]*$ ]]; then
        mem_mib="$SLURM_MEM_PER_NODE"
    elif [[ "${SLURM_MEM_PER_CPU:-}" =~ ^[1-9][0-9]*$ ]]; then
        mem_mib=$((SLURM_MEM_PER_CPU*cpus))
    elif [[ "$platform" == Darwin ]]; then
        mem_bytes="$(/usr/sbin/sysctl -n hw.memsize 2>/dev/null || true)"
        [[ "$mem_bytes" =~ ^[1-9][0-9]*$ ]] && mem_mib=$((mem_bytes/1024/1024)) || true
    elif [[ -r /proc/meminfo ]]; then
        detected="$(awk '$1 == "MemTotal:" {print int($2/1024); exit}' /proc/meminfo 2>/dev/null || true)"
        [[ "$detected" =~ ^[1-9][0-9]*$ ]] && mem_mib="$detected" || true
    fi
    [[ "$mem_mib" =~ ^[0-9]+$ ]] || mem_mib=0

    REFERENCE_BUILD_CPUS="$cpus"
    REFERENCE_BUILD_MEM_MIB="$mem_mib"
    BISMARK_INDEX_THREADS_PER_STRAND=$((cpus/2))
    (( BISMARK_INDEX_THREADS_PER_STRAND >= 1 )) || BISMARK_INDEX_THREADS_PER_STRAND=1
    (( BISMARK_INDEX_THREADS_PER_STRAND <= 8 )) || BISMARK_INDEX_THREADS_PER_STRAND=8
    BISMARK_COMBINED_INDEX_THREADS="$cpus"
    (( BISMARK_COMBINED_INDEX_THREADS >= 1 )) || BISMARK_COMBINED_INDEX_THREADS=1
    (( BISMARK_COMBINED_INDEX_THREADS <= 8 )) || BISMARK_COMBINED_INDEX_THREADS=8

    if (( cpus >= 4 && mem_mib >= 49152 )); then
        REFERENCE_PARALLEL_BISCUIT=1
        REFERENCE_PARALLEL_DOWNSTREAM=1
    else
        REFERENCE_PARALLEL_BISCUIT=0
        REFERENCE_PARALLEL_DOWNSTREAM=0
    fi

    echo "自适应 reference 资源：CPU=${REFERENCE_BUILD_CPUS}，内存=${REFERENCE_BUILD_MEM_MIB} MiB，Bismark CT/GA 每链线程=${BISMARK_INDEX_THREADS_PER_STRAND}，Combined 线程=${BISMARK_COMBINED_INDEX_THREADS}，双物种 BISCUIT 并行=${REFERENCE_PARALLEL_BISCUIT}，双物种 downstream 并行=$REFERENCE_PARALLEL_DOWNSTREAM"
}


setup_reference() {
    local species="$1" genome_source="$2" fasta="$3"
    local prefix assembled_ready bundle_ready

    echo
    echo "=== 标准 reference：$species ==="

    case "$species" in
        hg38)
            prefix=HG38
            assembled_ready="$HG38_ASSEMBLED_READY"
            bundle_ready="$HG38_BUNDLE_READY"
            ;;
        mm10)
            prefix=MM10
            assembled_ready="$MM10_ASSEMBLED_READY"
            bundle_ready="$MM10_BUNDLE_READY"
            ;;
        *)
            echo "ERROR: 不支持的受管 reference 物种： $species" >&2
            return 2
            ;;
    esac

    if (( bundle_ready == 1 )); then
        echo "$species reference bundle 已 READY，无需重建或刷新 manifest。"
        return 0
    fi

    if (( assembled_ready == 1 )); then
        echo "复用结构有效的 $species 组装 reference： $fasta"
    else
        echo "构建计划：组装后的 $species FASTA/FAI 未通过结构检查，重建 FASTA 与依赖 indexes。"
        assemble_reference "$species" "$genome_source" "$fasta" || return $?
        printf -v "${prefix}_BWA_READY" '%d' 0
        printf -v "${prefix}_BISCUIT_READY" '%d' 0
        printf -v "${prefix}_BISMARK_STANDARD_READY" '%d' 0
        printf -v "${prefix}_BISMARK_COMBINED_READY" '%d' 0
        printf -v "${prefix}_BISMARK_READY" '%d' 0
        printf -v "${prefix}_DOWNSTREAM_READY" '%d' 0
    fi

}

build_species_biscuit() {
    local species="$1" fasta prefix ready
    case "$species" in hg38) fasta="$HG38_REFERENCE"; prefix=HG38 ;; mm10) fasta="$MM10_REFERENCE"; prefix=MM10 ;; *) return 2 ;; esac
    eval "ready=\${${prefix}_BISCUIT_READY}"
    if (( ready == 1 )); then
        echo "复用完整 BISCUIT indexes： $fasta"
    else
        build_biscuit_index "$fasta"
    fi
}

run_biscuit_phase() {
    local hg_pid mm_pid hg_rc=0 mm_rc=0
    if (( REFERENCE_PARALLEL_BISCUIT == 1 && HG38_BISCUIT_READY == 0 && MM10_BISCUIT_READY == 0 )); then
        echo "并行构建 hg38 与 mm10 BISCUIT indexes。"
        build_species_biscuit hg38 & hg_pid=$!
        build_species_biscuit mm10 & mm_pid=$!
        REFERENCE_CHILD_PIDS="$hg_pid $mm_pid"
        if wait "$hg_pid"; then hg_rc=0; else hg_rc=$?; kill -TERM "$mm_pid" 2>/dev/null || true; fi
        if wait "$mm_pid"; then mm_rc=0; else mm_rc=$?; fi
        REFERENCE_CHILD_PIDS=""
        (( hg_rc == 0 && mm_rc == 0 )) || return 1
    else
        build_species_biscuit hg38 || return $?
        build_species_biscuit mm10 || return $?
    fi
}

build_species_bismark() {
    local species="$1" fasta prefix standard_ready combined_ready
    case "$species" in hg38) fasta="$HG38_REFERENCE"; prefix=HG38 ;; mm10) fasta="$MM10_REFERENCE"; prefix=MM10 ;; *) return 2 ;; esac
    eval "standard_ready=\${${prefix}_BISMARK_STANDARD_READY}"
    eval "combined_ready=\${${prefix}_BISMARK_COMBINED_READY}"

    if (( standard_ready == 1 )); then
        echo "复用标准 Bismark CT/GA indexes： $fasta.bismark"
    else
        build_bismark_standard_index "$fasta" || return $?
        standard_ready=1
        printf -v "${prefix}_BISMARK_STANDARD_READY" '%d' 1
    fi

    if (( combined_ready == 1 )); then
        echo "复用 Bismark Combined index： $fasta.bismark/Bisulfite_Genome/Combined"
    else
        build_bismark_combined_index "$fasta" || return $?
        combined_ready=1
        printf -v "${prefix}_BISMARK_COMBINED_READY" '%d' 1
    fi
    printf -v "${prefix}_BISMARK_READY" '%d' 1
}

run_bismark_phase() {
    build_species_bismark hg38 || return $?
    build_species_bismark mm10 || return $?
}

build_species_downstream() {
    local species="$1" fasta prefix ready
    case "$species" in hg38) fasta="$HG38_REFERENCE"; prefix=HG38 ;; mm10) fasta="$MM10_REFERENCE"; prefix=MM10 ;; *) return 2 ;; esac
    eval "ready=\${${prefix}_DOWNSTREAM_READY}"
    if (( ready == 1 )); then
        echo "复用 downstream TSS + single-CpG references： $(dirname "$fasta")"
    else
        build_downstream_references "$species" "$fasta"
    fi
}

run_downstream_phase() {
    local hg_pid mm_pid hg_rc=0 mm_rc=0
    if (( REFERENCE_PARALLEL_DOWNSTREAM == 1 && HG38_DOWNSTREAM_READY == 0 && MM10_DOWNSTREAM_READY == 0 )); then
        echo "并行构建 hg38 与 mm10 downstream references。"
        build_species_downstream hg38 & hg_pid=$!
        build_species_downstream mm10 & mm_pid=$!
        REFERENCE_CHILD_PIDS="$hg_pid $mm_pid"
        if wait "$hg_pid"; then hg_rc=0; else hg_rc=$?; kill -TERM "$mm_pid" 2>/dev/null || true; fi
        if wait "$mm_pid"; then mm_rc=0; else mm_rc=$?; fi
        REFERENCE_CHILD_PIDS=""
        (( hg_rc == 0 && mm_rc == 0 )) || return 1
    else
        build_species_downstream hg38 || return $?
        build_species_downstream mm10 || return $?
    fi
}

finalize_reference_species() {
    local species="$1" genome_source fasta prefix bundle_ready
    case "$species" in
        hg38) genome_source="$HG38_SOURCE"; fasta="$HG38_REFERENCE"; prefix=HG38 ;;
        mm10) genome_source="$MM10_SOURCE"; fasta="$MM10_REFERENCE"; prefix=MM10 ;;
        *) return 2 ;;
    esac
    eval "bundle_ready=\${${prefix}_BUNDLE_READY}"
    if (( bundle_ready == 0 )); then
        if ! validate_reference "$species" "$fasta"; then
            return 1
        fi
    fi
    echo "$species reference: READY"
}


slurm_resources_sufficient() {
    local cpus mem_mib

    [[ -n "${SLURM_JOB_ID:-}" ]] || {
        echo "未检测到活动的 Slurm allocation。"
        return 1
    }

    cpus="${SLURM_CPUS_ON_NODE:-}"
    if [[ ! "$cpus" =~ ^[1-9][0-9]*$ ]]; then
        cpus="${SLURM_CPUS_PER_TASK:-}"
    fi
    if [[ ! "$cpus" =~ ^[1-9][0-9]*$ ]]; then
        echo "Slurm job ${SLURM_JOB_ID}：无法确认已分配 CPU，将提交独立 reference 作业。"
        return 1
    fi

    if [[ "${SLURM_MEM_PER_NODE:-}" =~ ^[1-9][0-9]*$ ]]; then
        mem_mib="${SLURM_MEM_PER_NODE}"
    elif [[ "${SLURM_MEM_PER_CPU:-}" =~ ^[1-9][0-9]*$ ]]; then
        mem_mib=$(( SLURM_MEM_PER_CPU * cpus ))
    else
        echo "Slurm job ${SLURM_JOB_ID}：无法确认已分配内存，将提交独立 reference 作业。"
        return 1
    fi

    echo "当前 Slurm allocation：job=${SLURM_JOB_ID}，CPU=${cpus}，内存=${mem_mib} MiB"
    echo "Reference 资源下限：CPU>=${REFERENCE_MIN_CPUS}，内存>=${REFERENCE_MIN_MEM_MIB} MiB"

    if (( cpus >= REFERENCE_MIN_CPUS && mem_mib >= REFERENCE_MIN_MEM_MIB )); then
        echo "当前 Slurm allocation 资源充足，直接复用。"
        return 0
    fi

    echo "当前 Slurm allocation 低于 reference 资源下限，另行提交构建作业。"
    return 1
}

# 进度文本同时进 stdout（run_logged 路线收进 detail log）与控制 tty（stdout 被重定向时）。
progress_tee() {
    local data
    data="$(cat)"
    [[ -n "$data" ]] || return 0
    printf '%s\n' "$data"
    if [[ ! -t 1 && -w /dev/tty ]]; then
        printf '%s\n' "$data" 2>/dev/null > /dev/tty || true
    fi
    return 0
}

# 以当前用户的 sacct 记账终判；空记录重试，记账不可用时由调用方复验产物。
slurm_job_final_status() {
    local job_id="$1" retry_sleep="${2:-5}"
    local attempt row state exit_code
    command -v sacct >/dev/null 2>&1 || {
        echo "WARNING: sacct 不可用，无法读取作业 $job_id 的退出状态；由构建段终验确认结果。"
        return 0
    }
    state=""
    for attempt in 1 2 3 4 5 6; do
        row="$(sacct -j "$job_id" -X -n -P -o 'State,ExitCode' -u "$(id -un)" 2>/dev/null | head -n 1)" || row=""
        IFS='|' read -r state exit_code <<< "$row"
        state="${state//[[:space:]]/}"
        exit_code="${exit_code//[[:space:]]/}"
        [[ -n "$state" ]] && break
        state=""
        sleep "$retry_sleep"
    done
    if [[ -z "$state" ]]; then
        echo "WARNING: 作业 $job_id 尚无 accounting 记录；由构建段终验确认结果。"
        return 0
    fi
    echo "Slurm job $job_id 最终状态：state=$state exit=$exit_code"
    [[ "$state" == COMPLETED && "$exit_code" == 0:0 ]] && return 0
    return 1
}

# 轮询队列和日志后以 sacct 终判；Ctrl+C 只脱离等待，不取消作业。
wait_slurm_job_with_progress() {
    local job_id="$1" log_file="$2" label="$3" poll_seconds="${4:-30}"
    local waited=0 polls=0 last_state_key="" state_key queue_output
    local state runtime node reason last_log_lines=-1 log_lines saved_int_trap
    saved_int_trap="$(trap -p INT)"
    trap 'printf "job %s 仍在运行，已脱离等待；可 tail -f %s 观察或 scancel %s\n" "$job_id" "$log_file" "$job_id" >&2; exit 130' INT
    while :; do
        if ! queue_output="$(squeue -h -j "$job_id" -o '%T|%M|%N|%r' 2>&1)"; then
            [[ "$queue_output" == *"Invalid job id"* ]] && break
            printf '[%s] squeue 不可用，重试：%s\n' "$label" "$queue_output" | progress_tee
            sleep "$poll_seconds"
            waited=$((waited + poll_seconds))
            continue
        fi
        IFS='|' read -r state runtime node reason <<< "$queue_output"
        [[ -n "${state:-}" ]] || break
        [[ "${reason:-}" == "None" ]] && reason=""
        state_key="$state|$reason"
        if (( polls == 0 )) || [[ "$state_key" != "$last_state_key" ]] || (( polls % 10 == 0 )); then
            printf '[%s] waited %ss: state=%s runtime=%s node=%s%s' \
                "$label" "$waited" "$state" "$runtime" "$node" "${reason:+ reason=$reason}" | progress_tee
            last_state_key="$state_key"
        fi
        if [[ -f "$log_file" ]]; then
            log_lines="$(wc -l < "$log_file" 2>/dev/null | tr -d '[:space:]')"
            [[ "$log_lines" =~ ^[0-9]+$ ]] || log_lines=0
            if (( log_lines > last_log_lines )); then
                tail -n 2 -- "$log_file" 2>/dev/null | sed 's/^/    | /' | progress_tee
                last_log_lines="$log_lines"
            fi
        fi
        sleep "$poll_seconds"
        waited=$((waited + poll_seconds))
        polls=$((polls + 1))
    done
    if [[ -n "$saved_int_trap" ]]; then
        eval "$saved_int_trap"
    else
        trap - INT
    fi
    slurm_job_final_status "$job_id"
}

# 提交 reference 构建 SLURM 作业并轮询等待；失败时回读其日志尾部。
submit_slurm_reference_setup() {
    local sbatch_bin job_output job_id ref_log log_pattern

    sbatch_bin="$(command -v sbatch || true)"
    [[ -n "$sbatch_bin" && -x "$sbatch_bin" ]] || {
        echo "ERROR: 当前 Linux/HPC 环境找不到可执行的 sbatch。" >&2
        if [[ -n "${SLURM_JOB_ID:-}" ]]; then
            echo "当前位于 SLURM job/allocation（job=${SLURM_JOB_ID}），但该节点不能提交新的作业。" >&2
        fi
        echo "请回到可执行 sbatch 的 HPC 登录节点重新运行：" >&2
        echo "  bash core/doctor.sh reference" >&2
        return 127
    }

    mkdir -p "$REFERENCE_LOG_ROOT"
    log_pattern="$REFERENCE_LOG_ROOT/reference_setup_%j.log"

    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "当前 Slurm allocation 不满足或无法确认满足 reference 资源要求：job=${SLURM_JOB_ID}"
        echo "提交独立 reference 作业，当前作业等待其结束。"
    else
        echo "Linux/HPC 入口没有足够的 allocation，提交独立 reference 作业。"
    fi
    echo "Slurm 资源：${REFERENCE_SUBMIT_CPUS} CPU，${REFERENCE_SUBMIT_MEM} 内存"

    job_output="$($sbatch_bin \
        --parsable \
        --job-name=Alopex_reference_setup \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="$REFERENCE_SUBMIT_CPUS" \
        --mem="$REFERENCE_SUBMIT_MEM" \
        --chdir="$PIPELINE_ROOT" \
        --output="$log_pattern" \
        --error="$log_pattern" \
        --export=ALL \
        "$DOCTOR_SELF_PATH" reference)" || {
        echo "ERROR: reference 构建的 sbatch 提交失败。" >&2
        return 1
    }
    job_id="${job_output%%;*}"
    job_id="${job_id%%.*}"
    if ! [[ "$job_id" =~ ^[0-9]+$ ]]; then
        echo "ERROR: 无法从以下输出解析 Slurm job ID： $job_output" >&2
        return 1
    fi
    ref_log="$REFERENCE_LOG_ROOT/reference_setup_${job_id}.log"
    echo "SLURM job    : $job_id"
    echo "SLURM log    : $ref_log"

    if ! wait_slurm_job_with_progress "$job_id" "$ref_log" "reference build"; then
        echo "ERROR: Slurm reference 构建失败（job ${job_id}）。" >&2
        if [[ -f "$ref_log" ]]; then
            echo "--- Slurm 作业日志末尾： $ref_log ---" >&2
            tail -n 40 -- "$ref_log" >&2 || true
            echo "--- Slurm 作业日志结束 ---" >&2
        else
            echo "找不到 Slurm 作业日志；可列出候选： ls -t $REFERENCE_LOG_ROOT/reference_setup_*.log" >&2
        fi
        if [[ -n "${SLURM_JOB_ID:-}" ]]; then
            echo "如果集群禁止计算节点 nested sbatch，请回登录节点重新运行：" >&2
            echo "  bash core/doctor.sh reference" >&2
        fi
        return 1
    fi
    echo "Slurm reference 构建完成：job $job_id"
}

run_reference_build() {
    check_prerequisites_reference build || return $?

    prepare_reference_build_plan || return $?

    if (( HG38_ASSEMBLED_READY == 0 || MM10_ASSEMBLED_READY == 0 || HG38_DOWNSTREAM_READY == 0 || MM10_DOWNSTREAM_READY == 0 )); then
        check_sources_for_build_plan || return $?
    else
        echo "跳过 source FASTA/GTF 校验：本次仅需已有组装 reference/index。"
    fi

    detect_reference_build_resources || return $?

    setup_reference hg38 "$HG38_SOURCE" "$HG38_REFERENCE" || return $?
    setup_reference mm10 "$MM10_SOURCE" "$MM10_REFERENCE" || return $?
    run_biscuit_phase || return $?
    run_bismark_phase || return $?
    if (( HG38_BWA_READY == 0 )); then build_bwa_index "$HG38_REFERENCE" || return $?; fi
    if (( MM10_BWA_READY == 0 )); then build_bwa_index "$MM10_REFERENCE" || return $?; fi
    run_downstream_phase || return $?
    finalize_reference_species hg38 || return $?
    finalize_reference_species mm10 || return $?

    echo
    echo "SUCCESS: hg38 / mm10 FASTA、BISCUIT/Bismark/BWA index、TSS 与 single-CpG reference 均已物理验证。"
}


run_preflight_checks() {
    check_prerequisites_reference readonly

    if reference_bundle_ready hg38 "$HG38_REFERENCE" &&
       reference_bundle_ready mm10 "$MM10_REFERENCE"; then
        REFERENCES_READY=1
        echo "预检：hg38 / mm10 reference bundles 已 READY，复用已有 indexes。"
    else
        REFERENCES_READY=0
        echo "预检：至少一个 reference/index/downstream 组件需要修复。"
    fi
}


main_reference() {
    local platform

    run_preflight_checks

    if (( REFERENCES_READY == 1 )); then
        echo
        echo "SUCCESS: hg38 / mm10 reference bundles 已完整，无需构建。"
        return 0
    fi

    if (( CHECK_ONLY == 1 )); then
        echo "NOT READY: 至少一个标准 reference bundle 需要构建或修复。" >&2
        return 1
    fi

    platform="$(uname -s)"
    case "$platform" in
        Darwin)
            echo "Platform: macOS"
            with_transaction_lock_reference run_reference_build
            ;;
        Linux)
            echo "平台：Linux HPC"
            if slurm_resources_sufficient; then
                with_transaction_lock_reference run_reference_build
            else
                submit_slurm_reference_setup || return $?
                run_preflight_checks
                if (( REFERENCES_READY == 1 )); then
                    echo "SUCCESS: Slurm 构建后 reference bundles 已通过 READY 终验。"
                else
                    echo "ERROR: Slurm 构建后 reference bundles 仍未 READY。" >&2
                    return 1
                fi
            fi
            ;;
        *)
            echo "ERROR: 不支持的操作系统： $platform" >&2
            echo "支持的构建平台为 macOS（Darwin）和 Linux HPC。" >&2
            exit 2
            ;;
    esac
}


# reference 段入口：分派的唯一调用目标；flow 首行自带 set 严格模式。
run_reference_section() {
set -euo pipefail
CHECK_ONLY=0
if (( $# > 1 )); then
    echo "ERROR: core/doctor.sh reference 只接受 --check 或 --help。" >&2
    usage_reference >&2
    return 2
fi
if (( $# == 1 )); then
    case "$1" in
        --check)
            CHECK_ONLY=1
            ;;
        -h|--help|help)
            usage_reference
            return 0
            ;;
        *)
            echo "ERROR: core/doctor.sh reference 不接受参数：$1" >&2
            echo "hg38/mm10 由脚本统一管理；其它物种由用户自行构建。" >&2
            usage_reference >&2
            return 2
            ;;
    esac
fi
main_reference
}





DEMUX_ROOT="$PIPELINE_ROOT/core/demux_rs"
RUST_STATE_ROOT="$PIPELINE_ROOT/conda/rust"
RUST_LOG_ROOT="$RUST_STATE_ROOT/logs"
CARGO_HOME="$RUST_STATE_ROOT/cargo-home"

DEMUX_LOCK="$PIPELINE_ROOT/conda/.locks/demux-build.lock.d"

DEMUX_BIN="$RUST_STATE_ROOT/bin/demux_rs"

DEMUX_BUILD_CANDIDATE=""

cleanup_demux_build_candidate() {
    rm -f -- "$DEMUX_BIN.tmp.$$" 2>/dev/null || true
    DEMUX_BUILD_CANDIDATE=""
}

cleanup_demux_transaction() {
    cleanup_demux_build_candidate
    txn_release_lock "$DEMUX_LOCK"
}

demux_transaction_signal() {
    exit "$1"
}

with_transaction_lock_demux() {
    local rc
    txn_acquire_lock "$DEMUX_LOCK" "demux build transaction" || return $?

    trap 'cleanup_demux_transaction' EXIT
    trap 'demux_transaction_signal 130' INT
    trap 'demux_transaction_signal 143' TERM
    trap 'demux_transaction_signal 129' HUP

    if "$@"; then
        rc=0
    else
        rc=$?
    fi

    cleanup_demux_build_candidate
    trap - EXIT INT TERM HUP
    txn_release_lock "$DEMUX_LOCK"
    return "$rc"
}

SLURM_CPUS=8
SLURM_MEM="16G"

usage_demux() {
    cat <<'USAGE'
用法:
  bash core/doctor.sh demux [--check]

选项:
  --check
    只读校验已记录的 Rust 源码/二进制 identity 是否一致。

行为:
  macOS:
    本机拉取 Rust crates，然后在本机测试并构建 demux_rs。

  Linux HPC:
    先检查 Cargo.lock 依赖在 Pipeline 私有共享 Cargo 缓存中是否齐全。
    齐全时：已在 SLURM 作业内则立即离线构建，否则提交一个 8 CPU / 16 GiB 的 SLURM 构建作业。
    不齐全时：只允许在登录节点补齐缺失 crates；SLURM 计算节点内会失败，
    并提示用户回到登录节点重跑。
USAGE
}


CARGO_BIN="$CONTROL_ENV/bin/cargo"
RUSTC_BIN="$CONTROL_ENV/bin/rustc"

# 调用单模块 CLI 计算或校验 demux_rs 编译期源码身份。
run_demux_identity_cli() {
    env \
        -u PYTHONHOME -u PYTHONPATH -u PYTHONSTARTUP -u PYTHONUSERBASE -u VIRTUAL_ENV \
        PYTHONPATH="$PIPELINE_ROOT/core" \
        PYTHONNOUSERSITE=1 \
        "$PYTHON_BIN" -m dna_pipeline demux-identity "$@"
}

# 校验 control 环境与 Rust 工具链齐备、demux_rs 源码身份文件存在；build 模式下准备缓存目录。
check_prerequisites_demux() {
    [[ -x "$PYTHON_BIN" ]] || {
        echo "ERROR: 缺少 control 环境： $PYTHON_BIN" >&2
        echo "请先运行 core/doctor.sh；手动重建 Conda 可运行 core/doctor.sh conda。" >&2
        return 127
    }

    [[ -x "$CARGO_BIN" ]] || {
        echo "ERROR: 缺少 control Cargo： $CARGO_BIN" >&2
        echo "请先运行 core/doctor.sh；手动重建 Conda 可运行 core/doctor.sh conda。" >&2
        return 127
    }

    [[ -x "$RUSTC_BIN" ]] || {
        echo "ERROR: 缺少 control rustc： $RUSTC_BIN" >&2
        echo "请先运行 core/doctor.sh；手动重建 Conda 可运行 core/doctor.sh conda。" >&2
        return 127
    }

    [[ -f "$DEMUX_ROOT/Cargo.toml" && -f "$DEMUX_ROOT/Cargo.lock" ]] || {
        echo "ERROR: 以下目录缺少 Cargo.toml/Cargo.lock： $DEMUX_ROOT" >&2
        return 2
    }
    if [[ "${1:-build}" == "build" ]]; then
        mkdir -p "$RUST_STATE_ROOT" "$RUST_LOG_ROOT" "$CARGO_HOME"
    fi
}

fetch_dependencies() {
    echo "下载 Rust 依赖到 Pipeline 私有 Cargo 缓存： $CARGO_HOME"
    "$CARGO_BIN" fetch --manifest-path "$DEMUX_ROOT/Cargo.toml" --locked
}

check_offline_dependencies() {
    echo "检查共享 Cargo 缓存中的 Rust 依赖是否完整： $CARGO_HOME"
    echo "提示：网络文件系统上的校验可能数分钟没有输出。"
    local started="$SECONDS"
    "$CARGO_BIN" fetch \
        --manifest-path "$DEMUX_ROOT/Cargo.toml" \
        --locked \
        --offline || return $?
    echo "Rust 依赖缓存检查完成，用时 $((SECONDS - started))s。"
}

demux_source_revision() {
    run_demux_identity_cli revision "$DEMUX_ROOT"
}

verify_built_binary() {
    local expected_revision="$1" binary="${2:-$DEMUX_BIN}"
    [[ -x "$binary" ]] || {
        echo "ERROR: 构建后仍缺少 demux 二进制： $binary" >&2
        return 1
    }
    run_demux_identity_cli verify-binary "$binary" "$expected_revision" >/dev/null
}

# 在 demux 事务锁内复用构建缓存，路径按当前 RUST_STATE_ROOT 解析以保持沙箱隔离。
prepare_demux_build_candidate() {
    mkdir -p "$RUST_STATE_ROOT/candidates"
    DEMUX_BUILD_CANDIDATE="$RUST_STATE_ROOT/candidates/demux.target"
    export CARGO_TARGET_DIR="$DEMUX_BUILD_CANDIDATE"
}

publish_demux_binary() {
    local source_binary="$1" expected_revision="$2"
    local destination="$DEMUX_BIN" temporary="$DEMUX_BIN.tmp.$$"
    mkdir -p "$(dirname "$destination")"
    rm -f -- "$temporary"
    cp "$source_binary" "$temporary"
    chmod +x "$temporary"
    verify_built_binary "$expected_revision" "$temporary"
    mv -f -- "$temporary" "$destination"
}

# 在候选目录中离线/在线 cargo test + build demux_rs，校验 identity 后原子发布二进制。
build_demux() {
    local offline="$1"
    local jobs="${2:-}"
    local source_revision
    local -a cargo_args=(--release --locked)

    if [[ "$offline" == "1" ]]; then
        cargo_args+=(--offline)
    fi
    if [[ -n "$jobs" ]]; then
        cargo_args+=(-j "$jobs")
    fi

    cd "$DEMUX_ROOT"
    prepare_demux_build_candidate
    source_revision="$(demux_source_revision)"
    export DNA_PIPELINE_SOURCE_REVISION="$source_revision"

    echo "demux 源码身份： $DNA_PIPELINE_SOURCE_REVISION"
    echo "Rust 构建节点： $(hostname)"
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "SLURM job: $SLURM_JOB_ID"
    fi

    if [[ "$offline" == "1" ]]; then
        echo "正在离线测试 demux_rs……"
    else
        echo "正在本地测试 demux_rs……"
    fi
    if ! "$CARGO_BIN" test "${cargo_args[@]}"; then
        echo "ERROR: demux_rs 的 cargo test 失败，候选二进制不发布。" >&2
        return 1
    fi

    if [[ "$offline" == "1" ]]; then
        echo "正在离线构建 demux_rs……"
    else
        echo "正在本地构建 demux_rs……"
    fi
    if ! "$CARGO_BIN" build "${cargo_args[@]}"; then
        echo "ERROR: demux_rs 的 cargo build 失败，候选二进制不发布。" >&2
        return 1
    fi

    verify_built_binary "$source_revision" "$DEMUX_BUILD_CANDIDATE/release/demux_rs"
    publish_demux_binary "$DEMUX_BUILD_CANDIDATE/release/demux_rs" "$source_revision"
    cleanup_demux_build_candidate
    unset CARGO_TARGET_DIR
    printf 'Built: %s\n' "$DEMUX_BIN"
}

# 提交 8 CPU/16 GiB 的 SLURM 离线构建作业并轮询等待；成功后终验二进制 identity，失败时回读日志尾部。
submit_slurm_build() {
    local sbatch_bin job_output job_id job_log log_pattern
    sbatch_bin="$(command -v sbatch || true)"
    [[ -n "$sbatch_bin" && -x "$sbatch_bin" ]] || {
        echo "ERROR: 当前 Linux 主机找不到 sbatch。" >&2
        echo "请在具备 Slurm 客户端的 HPC login 节点运行 core/doctor.sh demux。" >&2
        return 127
    }

    log_pattern="$RUST_LOG_ROOT/demux_build_%j.log"
    echo "提交 Slurm 离线 demux 构建：8 CPU、16 GiB 内存"

    job_output="$($sbatch_bin \
        --parsable \
        --job-name=Alopex_demux_build \
        --nodes=1 \
        --ntasks=1 \
        --cpus-per-task="$SLURM_CPUS" \
        --mem="$SLURM_MEM" \
        --chdir="$PIPELINE_ROOT" \
        --output="$log_pattern" \
        --error="$log_pattern" \
        --export=ALL \
        "$DOCTOR_SELF_PATH" demux)" || {
        echo "ERROR: demux 构建的 sbatch 提交失败。" >&2
        return 1
    }
    job_id="${job_output%%;*}"
    job_id="${job_id%%.*}"
    if ! [[ "$job_id" =~ ^[0-9]+$ ]]; then
        echo "ERROR: 无法从以下输出解析 Slurm job ID： $job_output" >&2
        return 1
    fi
    job_log="$RUST_LOG_ROOT/demux_build_${job_id}.log"
    echo "SLURM job: $job_id"
    echo "SLURM log: $job_log"

    if ! wait_slurm_job_with_progress "$job_id" "$job_log" "demux build"; then
        echo "ERROR: Slurm demux 构建失败（job ${job_id}）。" >&2
        if [[ -f "$job_log" ]]; then
            echo "--- Slurm 作业日志末尾： $job_log ---" >&2
            tail -n 40 -- "$job_log" >&2 || true
            echo "--- Slurm 作业日志结束 ---" >&2
        else
            echo "找不到 Slurm 作业日志；可列出候选： ls -t $RUST_LOG_ROOT/demux_build_*.log" >&2
        fi
        return 1
    fi

    verify_built_binary "$(demux_source_revision)" || {
        echo "ERROR: Slurm 构建后 demux 二进制身份校验失败（job ${job_id}）。" >&2
        return 1
    }
    echo "Slurm demux 构建完成：job $job_id"
}

demux_fetch_and_build() {
    fetch_dependencies
    build_demux 0
}
demux_offline_build() {
    build_demux 1 "${SLURM_CPUS_PER_TASK:-$SLURM_CPUS}"
}
demux_fetch_only() {
    fetch_dependencies
}

# 主流程：--check 只读校验；否则按 macOS / Linux HPC 分路完成 fetch 与离线构建。
main_demux() {
    local platform revision
    platform="$(uname -s)"

    if (( CHECK_ONLY == 1 )); then
        check_prerequisites_demux check
        revision="$(demux_source_revision)"
        verify_built_binary "$revision"
        echo "READY: demux_rs 源码与二进制身份一致。"
        return 0
    fi

    case "$platform" in
        Darwin)
            echo "Platform: macOS"
            check_prerequisites_demux
            with_transaction_lock_demux demux_fetch_and_build
            ;;

        Linux)
            echo "平台：Linux HPC"
            check_prerequisites_demux

            if check_offline_dependencies; then
                echo "共享 Cargo 缓存中的 Rust 依赖完整。"

                if [[ -n "${SLURM_JOB_ID:-}" ]]; then
                    echo "检测到当前 Slurm 作业，在现有 allocation 内离线构建 demux_rs。"
                    with_transaction_lock_demux demux_offline_build
                else
                    echo "Rust 依赖已就绪，提交 Slurm 离线构建。"
                    submit_slurm_build
                fi
            else
                echo "共享 Cargo 缓存中的 Rust 依赖不完整： $CARGO_HOME" >&2

                if [[ -n "${SLURM_JOB_ID:-}" ]]; then
                    echo "ERROR: Slurm compute 作业内禁止下载缺失的 Rust 依赖。" >&2
                    echo "请在 HPC login 节点执行：" >&2
                    echo "  bash core/doctor.sh demux" >&2
                    exit 2
                fi

                echo "当前位于 HPC login 节点，下载缺失的 Rust 依赖。"
                with_transaction_lock_demux demux_fetch_only
                submit_slurm_build
            fi
            ;;

        *)
            echo "ERROR: 不支持的操作系统： $platform" >&2
            echo "支持的平台为 macOS（Darwin）和 Linux HPC。" >&2
            exit 2
            ;;
    esac
}


# 下载段入口：只预热 conda package_cache 与 cargo 依赖缓存；不创建/发布环境、不编译、不提交作业。
run_download_section() {
set -euo pipefail
CHECK_ONLY=0
export CARGO_HOME
if (( $# > 0 )); then
    echo "ERROR: core/doctor.sh download 不接受参数。" >&2
    doctor_usage >&2
    return 2
fi
local conda_bin conda_dir role specification warm_root
local -a specifier_args=()
require_specification "$ENVS_SPEC"
mkdir -p "$CONDA_ROOT" "$PACKAGE_CACHE" "$SETUP_HOME" "$SETUP_ENVS" "$SETUP_CACHE" "$SETUP_CONFIG" "$ENVS_BUILD_DIR"
warm_root="$CONDA_ROOT/.build/warm"
mkdir -p "$warm_root"
SEED_CONDA="$(resolve_seed_conda)"
SEED_CONDA_DIR="$(cd "$(dirname "$SEED_CONDA")" && pwd -P)"
if [[ -x "$CONDA_ROOT/current/control/bin/conda" ]]; then
    conda_bin="$CONDA_ROOT/current/control/bin/conda"
else
    conda_bin="$SEED_CONDA"
fi
conda_dir="$(cd "$(dirname "$conda_bin")" && pwd -P)"
run_conda_probe "$conda_bin" create --help 2>&1 | grep -q -- '--download-only' || {
    echo "ERROR: $conda_bin 不支持 conda create --download-only，无法预热包缓存。" >&2
    return 127
}
if run_conda_probe "$conda_bin" create --help 2>&1 | grep -q -- '--environment-specifier'; then
    specifier_args=(--environment-specifier environment-yaml)
fi
export PATH="$conda_dir:$HOST_PATH"
validate_seed_context
txn_acquire_lock "$ENVIRONMENT_LOCK" "conda cache warm" || return $?
trap cleanup_setup_transaction EXIT
for role in control biscuit bismark rastair notebook; do
    specification="$(extract_env_spec "$role")" || return $?
    echo "预热 Conda 包缓存： $role"
    env \
        "${CLEAN_ENVIRONMENT[@]}" \
        CONDA_ALWAYS_YES=true \
        CONDA_DEFAULT_THREADS=1 \
        CONDA_EXECUTE_THREADS=1 \
        CONDA_FETCH_THREADS=1 \
        CONDA_REMOTE_MAX_RETRIES=5 \
        CONDA_REPODATA_THREADS=1 \
        CONDA_VERIFY_THREADS=1 \
        PYTHONNOUSERSITE=1 \
        "$conda_bin" create \
            --download-only \
            "${specifier_args[@]}" \
            --prefix "$warm_root/$role" \
            --file "$specification"
done
txn_release_lock "$ENVIRONMENT_LOCK"
trap - EXIT
echo "预热 Rust 依赖缓存： cargo fetch"
with_transaction_lock_demux demux_fetch_only
echo "缓存预热完成： $PACKAGE_CACHE 与 ${CARGO_HOME}；后续 conda/demux 段可直接离线复用。"
}


# demux 段入口：分派的唯一调用目标；flow 首行自带 set 严格模式。
run_demux_section() {
set -euo pipefail
CHECK_ONLY=0
export CARGO_HOME
if (( $# > 1 )); then
    echo "ERROR: core/doctor.sh demux 只接受 --check 或 --help。" >&2
    usage_demux >&2
    return 2
fi
if (( $# == 1 )); then
    case "$1" in
        --check)
            CHECK_ONLY=1
            ;;
        -h|--help|help)
            usage_demux
            return 0
            ;;
        *)
            echo "ERROR: 未知参数： $1" >&2
            usage_demux >&2
            return 2
            ;;
    esac
fi
export RUSTC="$RUSTC_BIN"
main_demux
}



RUN_PIPELINE="$PIPELINE_ROOT/core/run_pipeline.sh"

BENCHMARK_ROOT="$PIPELINE_ROOT/resources/doctor_benchmark"
BENCHMARK_JSON="$BENCHMARK_ROOT/benchmark.json"

DOCTOR_TMP_PARENT="$PIPELINE_ROOT/logs"

DOCTOR_SANDBOX_RETENTION_DAYS=14


# 清理 Doctor 初始化期间被中断的本轮沙箱。
doctor_early_cleanup() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ -n "${TMP_ROOT:-}" && "$TMP_ROOT" == "$DOCTOR_TMP_PARENT/doctor."* && -e "$TMP_ROOT" ]]; then
        find "$TMP_ROOT" -type d -exec chmod u+rwx {} + 2>> "$DETAIL_LOG" || true
        find "$TMP_ROOT" -type f -exec chmod u+rw {} + 2>> "$DETAIL_LOG" || true
        rm -rf -- "$TMP_ROOT" >> "$DETAIL_LOG" 2>&1 || true
    fi
    exit "$status"
}

EXECUTOR=""
CURRENT_STEP="initialization"
FINAL_STATUS="NOT READY"
DOCTOR_TOTAL_STEPS=9
DOCTOR_STEP=0
DOCTOR_STEP_STARTED=0
DOCTOR_ACTIVE_COMMAND_PID=""
DOCTOR_ACTIVE_ROUTE_PIDS=""
DOCTOR_ACTIVE_ROUTE_JOB_IDS=""


# 有界终止的进度与告警统一进 detail 日志。
doctor_terminate_pid_bounded() {
    txn_terminate_pid_bounded "$1" "$2" "${DOCTOR_TERMINATION_GRACE_SECONDS:-15}" "${3:-0}" \
        >> "$DETAIL_LOG" 2>&1
}
doctor_signal_process_tree() {
    txn_signal_process_tree "$1" "$2"
}
stop_active_command() {
    local active_pid="${DOCTOR_ACTIVE_COMMAND_PID:-}"
    DOCTOR_ACTIVE_COMMAND_PID=""
    [[ "$active_pid" =~ ^[0-9]+$ ]] || return 0
    doctor_terminate_pid_bounded "$active_pid" "active Doctor command"
}
stop_active_routes() {
    local active_pid
    local job_id
    for job_id in ${DOCTOR_ACTIVE_ROUTE_JOB_IDS:-}; do
        [[ "$job_id" =~ ^[0-9]+$ ]] || continue
        scancel "$job_id" >> "$DETAIL_LOG" 2>&1 || true
    done
    DOCTOR_ACTIVE_ROUTE_JOB_IDS=""
    for active_pid in ${DOCTOR_ACTIVE_ROUTE_PIDS:-}; do
        [[ "$active_pid" =~ ^[0-9]+$ ]] || continue
        if txn_pid_is_running "$active_pid"; then
            printf '向活动 Doctor route（pid=%s）发送 TERM\n' "$active_pid" >> "$DETAIL_LOG"
            doctor_signal_process_tree "$active_pid" TERM
        fi
    done
    for active_pid in ${DOCTOR_ACTIVE_ROUTE_PIDS:-}; do
        [[ "$active_pid" =~ ^[0-9]+$ ]] || continue
        doctor_terminate_pid_bounded "$active_pid" "active Doctor route" 1
    done
    DOCTOR_ACTIVE_ROUTE_PIDS=""
}
doctor_refresh_active_route_pids() {
    DOCTOR_ACTIVE_ROUTE_PIDS=""
    [[ -z "${ROUTE_A_RC:-}" && "${ROUTE_A_PID:-}" =~ ^[0-9]+$ ]] && DOCTOR_ACTIVE_ROUTE_PIDS="$ROUTE_A_PID"
    if [[ -z "${ROUTE_B_RC:-}" && "${ROUTE_B_PID:-}" =~ ^[0-9]+$ ]]; then
        DOCTOR_ACTIVE_ROUTE_PIDS="${DOCTOR_ACTIVE_ROUTE_PIDS:+$DOCTOR_ACTIVE_ROUTE_PIDS }$ROUTE_B_PID"
    fi
}
doctor_fail_fast_routes() {
    local now
    now="$(date +%s)"
    if [[ -n "${ROUTE_A_RC:-}" ]] && (( ROUTE_A_RC != 0 )) && [[ -z "${ROUTE_B_RC:-}" ]]; then
        printf 'Route A 失败，rc=%s；请求取消 Route B。\n' "$ROUTE_A_RC" >> "$DETAIL_LOG"
        [[ "${ROUTE_B_JOB_ID:-}" =~ ^[0-9]+$ ]] && scancel "$ROUTE_B_JOB_ID" >> "$DETAIL_LOG" 2>&1 || true
        doctor_terminate_pid_bounded "$ROUTE_B_PID" "Route B after Route A failure"
        ROUTE_B_RC=143
        ROUTE_B_FINISHED="$now"
        ROUTE_B_CANCELLED=1
    elif [[ -n "${ROUTE_B_RC:-}" ]] && (( ROUTE_B_RC != 0 )) && [[ -z "${ROUTE_A_RC:-}" ]]; then
        printf 'Route B 失败，rc=%s；请求取消 Route A。\n' "$ROUTE_B_RC" >> "$DETAIL_LOG"
        [[ "${ROUTE_A_JOB_ID:-}" =~ ^[0-9]+$ ]] && scancel "$ROUTE_A_JOB_ID" >> "$DETAIL_LOG" 2>&1 || true
        doctor_terminate_pid_bounded "$ROUTE_A_PID" "Route A after Route B failure"
        ROUTE_A_RC=143
        ROUTE_A_FINISHED="$now"
        ROUTE_A_CANCELLED=1
    fi
    doctor_refresh_active_route_pids
}

# Doctor 小型路线在 allocation 内本地运行；login 节点按路线提交独立的小型整路线作业。
start_doctor_route() {
    local route="$1" project="$2" route_log="$3" job_var="$4" pid_var="$5"
    local worker job_output job_id worker_mem
    case "$route" in
        hg38-cabernet-bismark) worker_mem=16G ;;
        mm10-srd-biscuit|hg38-taps-rastair) worker_mem=24G ;;
        *) return 2 ;;
    esac
    if [[ "$SYSTEM_NAME" == Linux && -z "${SLURM_JOB_ID:-}" ]]; then
        worker="$TMP_ROOT/${route}.worker.sh"
        cat > "$worker" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd $(printf '%q' "$project")
exec env DNA_PIPELINE_DATA_MODE=test DNA_PIPELINE_TEST_ROUTE=$(printf '%q' "$route") DNA_PIPELINE_DOCTOR_LOCAL_CORES=2 bash $(printf '%q' "$RUN_PIPELINE")
EOF
        chmod +x "$worker"
        job_output="$(sbatch --parsable --job-name="Alopex_doctor_${route}" --nodes=1 --ntasks=1 --cpus-per-task=2 --mem="$worker_mem" --time=00:20:00 --chdir="$PIPELINE_ROOT" --output="$route_log" --error="$route_log" "$worker")" || return $?
        job_id="${job_output%%;*}"; job_id="${job_id%%.*}"
        [[ "$job_id" =~ ^[0-9]+$ ]] || return 1
        printf -v "$job_var" '%s' "$job_id"
        wait_slurm_job_with_progress "$job_id" "$route_log" "$route" >> "$DETAIL_LOG" 2>&1 &
    else
        printf -v "$job_var" '%s' ""
        env DNA_PIPELINE_DATA_MODE=test DNA_PIPELINE_TEST_ROUTE="$route" DNA_PIPELINE_DOCTOR_LOCAL_CORES=2 \
            bash -c 'cd "$1" && exec bash "$2"' _ "$project" "$RUN_PIPELINE" > "$route_log" 2>&1 &
    fi
    printf -v "$pid_var" '%s' "$!"
}
remove_doctor_sandbox() {
    local expected_prefix="$DOCTOR_TMP_PARENT/doctor."
    [[ -n "${TMP_ROOT:-}" && "$TMP_ROOT" == "$expected_prefix"* && "$TMP_ROOT" != "$DOCTOR_TMP_PARENT" ]] || {
        printf '拒绝清理不安全的 Doctor 沙箱路径： %s\n' "${TMP_ROOT:-<empty>}" >> "$DETAIL_LOG"
        return 2
    }
    [[ -e "$TMP_ROOT" ]] || return 0

    find "$TMP_ROOT" -type d -exec chmod u+rwx {} + 2>> "$DETAIL_LOG" || true
    find "$TMP_ROOT" -type f -exec chmod u+rw {} + 2>> "$DETAIL_LOG" || true
    rm -rf -- "$TMP_ROOT" >> "$DETAIL_LOG" 2>&1
    [[ ! -e "$TMP_ROOT" ]]
}
summary() {
    printf '%s\n' "$*"
    printf '%s\n' "$*" >> "$FINAL_LOG"
}
detail_section() { printf '\n===== %s =====\n' "$1" >> "$DETAIL_LOG"; }

# 只为完整阶段编号；构建修复与 manifest 刷新作为阶段内部命令记录。
doctor_step_start() {
    DOCTOR_STEP=$((DOCTOR_STEP + 1))
    CURRENT_STEP="$1"
    DOCTOR_STEP_STARTED="$(date +%s)"
    detail_section "$CURRENT_STEP"
    summary "Step $DOCTOR_STEP/$DOCTOR_TOTAL_STEPS START: $CURRENT_STEP"
}
doctor_step_ok() {
    local note="${1:+; $1}"
    summary "Step $DOCTOR_STEP/$DOCTOR_TOTAL_STEPS OK: $CURRENT_STEP ($(( $(date +%s) - DOCTOR_STEP_STARTED ))s$note)"
}

# 内部命令的原生输出与退出码仅进 detail；保留活动 PID 以便中断时收尾。
run_logged() {
    local label="$1" rc started; shift
    started="$(date +%s)"; detail_section "$label"
    { printf '+ '; printf '%q ' "$@"; printf '\n'; } >> "$DETAIL_LOG"
    "$@" >> "$DETAIL_LOG" 2>&1 & DOCTOR_ACTIVE_COMMAND_PID=$!
    wait "$DOCTOR_ACTIVE_COMMAND_PID"; rc=$?
    DOCTOR_ACTIVE_COMMAND_PID=""
    printf '命令结束：%s（rc=%s；%ss）\n' "$label" "$rc" "$(( $(date +%s) - started ))" >> "$DETAIL_LOG"
    return "$rc"
}
# Doctor 全局收尾器：停活动子进程、按成败保留/清理沙箱并输出最终状态。
finish_doctor() {
    local status=$? cleanup_status="已删除"
    trap - EXIT INT TERM HUP
    stop_active_command
    stop_active_routes
    if (( status != 0 )); then
        summary "Step $DOCTOR_STEP/$DOCTOR_TOTAL_STEPS FAIL: $CURRENT_STEP （rc=${status}；详见 detail 日志）"
    fi
    if (( status == 0 )); then
        if ! remove_doctor_sandbox; then
            cleanup_status="删除失败（保留：${TMP_ROOT}）"
            summary "[FAIL] Doctor 无法删除已完成的沙箱： $TMP_ROOT"
            status=1
        fi
    else
        cleanup_status="为失败诊断保留：$TMP_ROOT"
    fi
    if (( status == 0 )); then FINAL_STATUS=READY; else FINAL_STATUS="NOT READY"; fi
    summary ""
    summary "FINAL STATUS: $FINAL_STATUS"
    summary "临时沙箱： $cleanup_status"
    exit "$status"
}
ensure_builder_ready() {
    local label="$1"
    local section="$2"

    doctor_step_start "$label"
    if run_logged "$label check" bash "$DOCTOR_SELF_PATH" "$section" --check; then
        doctor_step_ok
        return 0
    fi

    summary "[REPAIR] $label 尚未 READY，执行 $section 构建"
    run_logged "$label repair" bash "$DOCTOR_SELF_PATH" "$section" || {
        summary "[FAIL] $label 修复失败"
        return 1
    }
    run_logged "$label re-check" bash "$DOCTOR_SELF_PATH" "$section" --check || {
        summary "[FAIL] $label 修复完成，但 READY 复验失败"
        return 1
    }
    doctor_step_ok "已修复并通过复验"
}

CONTROL_PYTHON="$CONTROL_ENV/bin/python"
BISMARK_SAMTOOLS="$BISMARK_ENV/bin/samtools"
BISCUIT_SAMTOOLS="$BISCUIT_ENV/bin/samtools"

# 在临时沙箱内摆放一条 smoke 项目的目录骨架、FASTQ、Barcode Map 与路线 config。
stage_smoke_project() {
    local project="$1" route="$2" reference="$3"
    PYTHONPATH="$PIPELINE_ROOT/core" "$CONTROL_PYTHON" -m dna_pipeline benchmark-stage \
        --project "$project" --reference "$reference" --route "$route" || return 1
}

run_manifest_refresh() {
    local project="$1"
    (
        unset PYTHONHOME
        export PYTHONNOUSERSITE=1
        export PYTHONPATH="$PIPELINE_ROOT/core"
        exec "$CONTROL_PYTHON" -m dna_pipeline refresh-manifest --project "$project"
    )
}

# 对照 benchmark 基准逐项校验路线的 sealed 输出与 demux 不变量。
verify_route_outputs() {
    local project="$1" route="$2"
    PYTHONPATH="$PIPELINE_ROOT/core" "$CONTROL_PYTHON" -m dna_pipeline benchmark-verify \
        --project "$project" --manifest "$BENCHMARK_JSON" --route "$route" || return 1
    cp "$project/04_logs/doctor_verification.json" "$DOCTOR_LOG_ROOT/doctor.$DOCTOR_RUN_ID.$route.json"
}

# doctor 段入口：分派的唯一调用目标；flow 首行自带 set 严格模式。
run_doctor_section() {
set -uo pipefail
DOCTOR_RUN_ID="$(date '+%Y%m%dT%H%M%S').$$"
DOCTOR_LOG_ROOT="${DNA_PIPELINE_DOCTOR_LOG_DIR:-$PIPELINE_ROOT/logs/doctor}"
mkdir -p "$DOCTOR_LOG_ROOT" "$DOCTOR_TMP_PARENT"
find "$DOCTOR_TMP_PARENT" -mindepth 1 -maxdepth 1 -type d -name 'doctor.*' \
    -mtime "+$DOCTOR_SANDBOX_RETENTION_DAYS" \
    -exec rm -rf -- {} + 2>/dev/null || true
TMP_ROOT="$(mktemp -d "$DOCTOR_TMP_PARENT/doctor.${DOCTOR_RUN_ID}.XXXXXX")" || exit 1
DETAIL_LOG="$DOCTOR_LOG_ROOT/doctor.${DOCTOR_RUN_ID}.detail.log"
FINAL_LOG="$DOCTOR_LOG_ROOT/doctor.${DOCTOR_RUN_ID}.summary.log"
LATEST_DETAIL_LOG="$DOCTOR_LOG_ROOT/latest.detail.log"
LATEST_SUMMARY_LOG="$DOCTOR_LOG_ROOT/latest.summary.log"
SMOKE_HG38="$TMP_ROOT/smoke_hg38_cabernet_bismark"
SMOKE_MM10="$TMP_ROOT/smoke_mm10_srd_biscuit"
: > "$DETAIL_LOG"; : > "$FINAL_LOG"
ln -sfn "$(basename "$DETAIL_LOG")" "$LATEST_DETAIL_LOG"
ln -sfn "$(basename "$FINAL_LOG")" "$LATEST_SUMMARY_LOG"
summary "Alopex Doctor"
summary "Pipeline : $PIPELINE_ROOT"
summary "Detail log : $DETAIL_LOG"
summary "Summary log: $FINAL_LOG"
trap doctor_early_cleanup EXIT
trap 'exit 130' INT TERM HUP
SYSTEM_NAME="$(uname -s)"
MACHINE_NAME="$(uname -m)"
HOST_NAME="$(hostname 2>/dev/null || uname -n)"
DOCTOR_TERMINATION_GRACE_SECONDS="${DNA_PIPELINE_TERMINATION_GRACE_SECONDS:-15}"
[[ "$DOCTOR_TERMINATION_GRACE_SECONDS" =~ ^[1-9][0-9]*$ ]] || DOCTOR_TERMINATION_GRACE_SECONDS=15
trap finish_doctor EXIT
printf 'Host: %s\nPlatform: %s/%s\n' "$HOST_NAME" "$SYSTEM_NAME" "$MACHINE_NAME" >> "$DETAIL_LOG"
summary "Route A  : hg38 / Cabernet / Bismark / 真实 hg38_cabernet raw + 合成对照"
summary "Route B  : mm10 / SRD / BISCUIT / 真实 mm10_srd Cabernet DNA + 合成 RNA 布局/对照"
summary "Route C  : hg38 / Cabernet-TAPS+ / Rastair / 独立合成真值"
summary ""
doctor_step_start "平台 / executor 检测"
case "$SYSTEM_NAME" in
    Darwin)
        EXECUTOR="local"
        ;;
    Linux)
        if ! command -v sbatch >/dev/null 2>&1 || ! command -v squeue >/dev/null 2>&1 || ! command -v scancel >/dev/null 2>&1; then
            summary "[FAIL] Linux Doctor 需要 HPC Slurm 命令（sbatch + squeue + scancel）。"
            return 1
        fi
        EXECUTOR="slurm"
        ;;
    *)
        summary "[FAIL] 不支持的操作系统： $SYSTEM_NAME"
        return 1
        ;;
esac
doctor_step_ok "$SYSTEM_NAME/${MACHINE_NAME}；生产 executor=$EXECUTOR"
ensure_builder_ready "Conda release" conda || exit 1
ensure_builder_ready "Rust demux_rs" demux || exit 1
ensure_builder_ready "hg38/mm10 reference bundle" reference || exit 1
doctor_step_start "benchmark 输入 / raw FASTQ 完整性"
for required in "$CONTROL_PYTHON" "$BISMARK_SAMTOOLS" "$BISCUIT_SAMTOOLS" "$RUN_PIPELINE"; do
    [[ -e "$required" ]] || {
        summary "[FAIL] 修复后仍缺少必需运行组件： $required"
        return 1
    }
done
for reference in "$HG38_REFERENCE" "$MM10_REFERENCE"; do
    [[ -s "$reference" && -s "$reference.fai" ]] || {
        summary "[FAIL] Doctor smoke reference 尚未 READY： $reference"
        return 1
    }
done
[[ -d "$BENCHMARK_ROOT" && -s "$BENCHMARK_JSON" ]] || {
    summary "[FAIL] 缺少固定 Doctor benchmark： $BENCHMARK_ROOT"
    return 1
}
PYTHONPATH="$PIPELINE_ROOT/core" "$CONTROL_PYTHON" -m dna_pipeline benchmark-validate \
    --root "$BENCHMARK_ROOT" --manifest "$BENCHMARK_JSON" >> "$DETAIL_LOG" 2>&1 || {
    summary "[FAIL] 固定 Doctor benchmark 完整性校验失败"
    return 1
}
doctor_step_ok "2100 对真实 reads / 4 个 barcode；SHA256 与配对校验通过"
doctor_step_start "Route A/B 项目准备与 manifest 刷新"
stage_smoke_project "$SMOKE_HG38" hg38-cabernet-bismark "$HG38_REFERENCE" >> "$DETAIL_LOG" 2>&1 || {
    summary "[FAIL] 无法准备 hg38/Cabernet/Bismark smoke 项目"
    return 1
}
stage_smoke_project "$SMOKE_MM10" mm10-srd-biscuit "$MM10_REFERENCE" >> "$DETAIL_LOG" 2>&1 || {
    summary "[FAIL] 无法准备 mm10/SRD/BISCUIT smoke 项目"
    return 1
}
run_logged "Route A manifest refresh" run_manifest_refresh "$SMOKE_HG38" || {
    summary "[FAIL] hg38/Cabernet/Bismark manifest 刷新失败"
    return 1
}
run_logged "Route B manifest refresh" run_manifest_refresh "$SMOKE_MM10" || {
    summary "[FAIL] mm10/SRD/BISCUIT manifest 刷新失败"
    return 1
}
doctor_step_ok
doctor_step_start "Route A/B 并行 workflow"
ROUTE_A_LOG="$TMP_ROOT/route_a.hg38-cabernet-bismark.log"
ROUTE_B_LOG="$TMP_ROOT/route_b.mm10-srd-biscuit.log"
HG38_START_EPOCH="$(date +%s)"
MM10_START_EPOCH="$HG38_START_EPOCH"
printf '+ route A: hg38/Cabernet/Bismark\n' >> "$DETAIL_LOG"
ROUTE_A_JOB_ID=""
start_doctor_route hg38-cabernet-bismark "$SMOKE_HG38" "$ROUTE_A_LOG" ROUTE_A_JOB_ID ROUTE_A_PID || {
    summary "[FAIL] 无法启动 Route A"
    return 1
}
printf '+ route B: mm10/SRD/BISCUIT\n' >> "$DETAIL_LOG"
ROUTE_B_JOB_ID=""
start_doctor_route mm10-srd-biscuit "$SMOKE_MM10" "$ROUTE_B_LOG" ROUTE_B_JOB_ID ROUTE_B_PID || {
    [[ "$ROUTE_A_JOB_ID" =~ ^[0-9]+$ ]] && scancel "$ROUTE_A_JOB_ID" >> "$DETAIL_LOG" 2>&1 || true
    doctor_terminate_pid_bounded "$ROUTE_A_PID" "Route A after Route B launch failure"
    summary "[FAIL] 无法启动 Route B"
    return 1
}
DOCTOR_ACTIVE_ROUTE_PIDS="$ROUTE_A_PID $ROUTE_B_PID"
DOCTOR_ACTIVE_ROUTE_JOB_IDS="${ROUTE_A_JOB_ID:+$ROUTE_A_JOB_ID }${ROUTE_B_JOB_ID:-}"
if [[ "$SYSTEM_NAME" == Linux && -z "${SLURM_JOB_ID:-}" ]]; then
    DOCTOR_ROUTE_EXECUTION_LABEL="two Slurm workers (2 CPU each; Bismark 16 GiB / BISCUIT 24 GiB; local Snakemake with mem-budgeted scheduling)"
elif [[ "$SYSTEM_NAME" == Linux ]]; then
    DOCTOR_ROUTE_EXECUTION_LABEL="current allocation (2 local cores per route)"
else
    DOCTOR_ROUTE_EXECUTION_LABEL="local (2 cores per route)"
fi
ROUTE_A_RC=""; ROUTE_B_RC=""; ROUTE_A_FINISHED=""; ROUTE_B_FINISHED=""
ROUTE_A_CANCELLED=0; ROUTE_B_CANCELLED=0
while [[ -z "$ROUTE_A_RC" || -z "$ROUTE_B_RC" ]]; do
    ROUTE_NOW="$(date +%s)"
    if [[ -z "$ROUTE_A_RC" ]] && ! kill -0 "$ROUTE_A_PID" 2>/dev/null; then
        wait "$ROUTE_A_PID"; ROUTE_A_RC=$?; ROUTE_A_FINISHED="$ROUTE_NOW"
        printf 'Route A（hg38/Cabernet/Bismark）结束，rc=%s\n' "$ROUTE_A_RC" >> "$DETAIL_LOG"
    fi
    if [[ -z "$ROUTE_B_RC" ]] && ! kill -0 "$ROUTE_B_PID" 2>/dev/null; then
        wait "$ROUTE_B_PID"; ROUTE_B_RC=$?; ROUTE_B_FINISHED="$ROUTE_NOW"
        printf 'Route B（mm10/SRD/BISCUIT）结束，rc=%s\n' "$ROUTE_B_RC" >> "$DETAIL_LOG"
    fi
    doctor_fail_fast_routes
    [[ -n "$ROUTE_A_RC" && -n "$ROUTE_B_RC" ]] || sleep 2
done
DOCTOR_ACTIVE_ROUTE_JOB_IDS=""
HG38_WALL_SECONDS=$((ROUTE_A_FINISHED-HG38_START_EPOCH))
MM10_WALL_SECONDS=$((ROUTE_B_FINISHED-MM10_START_EPOCH))
{
    printf '\n===== Route A: hg38/Cabernet/Bismark =====\n'
    cat "$ROUTE_A_LOG"
    printf '\n===== Route B: mm10/SRD/BISCUIT =====\n'
    cat "$ROUTE_B_LOG"
} >> "$DETAIL_LOG"
if (( ROUTE_A_RC != 0 || ROUTE_B_RC != 0 )); then
    append_controller_failure_tail() {
        local project="$1" route_label="$2" controller_log found=0
        for controller_log in "$project"/04_logs/controller/*.log; do
            [[ -f "$controller_log" ]] || continue
            found=1
            {
                printf '\n--- %s controller 失败日志末尾：%s ---\n' "$route_label" "$controller_log"
                tail -n 160 "$controller_log" 2>/dev/null || true
            } >> "$DETAIL_LOG"
        done
        (( found == 1 )) || printf '\n--- %s controller 未生成失败日志 ---\n' "$route_label" >> "$DETAIL_LOG"
    }
    (( ROUTE_A_RC == 0 )) || append_controller_failure_tail "$SMOKE_HG38" "Route A"
    (( ROUTE_B_RC == 0 )) || append_controller_failure_tail "$SMOKE_MM10" "Route B"
    if (( ROUTE_A_CANCELLED == 1 )); then
        summary "[CANCELLED] Route B 失败后已请求取消 Route A"
    elif (( ROUTE_A_RC != 0 )); then
        summary "[FAIL] Route A hg38/Cabernet/Bismark workflow 失败（rc=${ROUTE_A_RC}）"
    fi
    if (( ROUTE_B_CANCELLED == 1 )); then
        summary "[CANCELLED] Route A 失败后已请求取消 Route B"
    elif (( ROUTE_B_RC != 0 )); then
        summary "[FAIL] Route B mm10/SRD/BISCUIT workflow 失败（rc=${ROUTE_B_RC}）"
    fi
    return 1
fi
doctor_step_ok "A=${HG38_WALL_SECONDS}s; B=${MM10_WALL_SECONDS}s; $DOCTOR_ROUTE_EXECUTION_LABEL"
doctor_step_start "Route A/B 输出 / mapping 校验"
verify_route_outputs "$SMOKE_HG38" hg38-cabernet-bismark >> "$DETAIL_LOG" 2>&1 || {
    summary "[FAIL] Route A 输出/demux 校验失败"
    return 1
}
verify_route_outputs "$SMOKE_MM10" mm10-srd-biscuit >> "$DETAIL_LOG" 2>&1 || {
    summary "[FAIL] Route B 输出/demux 校验失败"
    return 1
}
doctor_step_ok
doctor_step_start "Route C TAPS / Rastair 项目准备、workflow 与合成真值校验"
SMOKE_TAPS="$TMP_ROOT/taps"
ROUTE_TAPS_LOG="$TMP_ROOT/route_c.hg38-taps-rastair.log"
stage_smoke_project "$SMOKE_TAPS" hg38-taps-rastair "$HG38_REFERENCE" >> "$DETAIL_LOG" 2>&1 || { summary "[FAIL] 无法准备 Route C"; return 1; }
run_logged "Route C manifest refresh" run_manifest_refresh "$SMOKE_TAPS" || { summary "[FAIL] Route C manifest 刷新失败"; return 1; }
start_doctor_route hg38-taps-rastair "$SMOKE_TAPS" "$ROUTE_TAPS_LOG" ROUTE_TAPS_JOB_ID ROUTE_TAPS_PID || { summary "[FAIL] 无法启动 Route C"; return 1; }
DOCTOR_ACTIVE_ROUTE_PIDS="$ROUTE_TAPS_PID"
DOCTOR_ACTIVE_ROUTE_JOB_IDS="${ROUTE_TAPS_JOB_ID:-}"
wait "$ROUTE_TAPS_PID"; ROUTE_TAPS_RC=$?
DOCTOR_ACTIVE_ROUTE_PIDS=""
DOCTOR_ACTIVE_ROUTE_JOB_IDS=""
cat "$ROUTE_TAPS_LOG" >> "$DETAIL_LOG"
(( ROUTE_TAPS_RC == 0 )) || { summary "[FAIL] TAPS workflow 失败： $ROUTE_TAPS_LOG"; return 1; }
verify_route_outputs "$SMOKE_TAPS" hg38-taps-rastair >> "$DETAIL_LOG" 2>&1 || { summary "[FAIL] TAPS 真值校验失败"; return 1; }
doctor_step_ok
CURRENT_STEP="complete"
return 0
}

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    return 0
fi

SECTION="${1:-doctor}"
case "$SECTION" in
    conda|reference|demux|download)
        shift
        ;;
    doctor)
        ;;
    -h|--help)
        doctor_usage
        exit 0
        ;;
    *)
        echo "ERROR: core/doctor.sh 只接受 conda|reference|demux|download 段或无参数（完整 doctor）。" >&2
        doctor_usage >&2
        exit 2
        ;;
esac

case "$SECTION" in
    conda) run_conda_section "$@" ;;
    reference) run_reference_section "$@" ;;
    demux) run_demux_section "$@" ;;
    download) run_download_section "$@" ;;
    doctor) run_doctor_section ;;
esac
exit $?

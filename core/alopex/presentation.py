"""仅观测本次运行的输出、进程和队列；不参与科学计算、调度或完成判定。"""
from __future__ import annotations

from collections import Counter, deque
from contextlib import nullcontext
from dataclasses import dataclass, field
from datetime import datetime
import html
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import threading

_LOADED_MODULE_CODE = sys._getframe().f_code

import time
from typing import TextIO


ANSI = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
CHILD_ENV = "DNA_PIPELINE_PRESENTATION_CHILD"
SNAKEMAKE_JOB_ID = r"(?:[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}|\d+)"


def observe_slurm_group_completion(workflow) -> None:
    """补记已完成 Slurm 分组的 UUID，和插件的分组提交日志对应。"""
    if getattr(getattr(workflow, 'executor_plugin', None), 'cli_prefix', None) != 'slurm':
        return
    executor = getattr(workflow.scheduler, '_executor', None)
    completed = getattr(executor, 'handle_job_success', None)
    if completed is None:
        return

    def report(job):
        result = completed(job)
        try:
            if job.is_group():
                print(f'[Alopex group finished] {job.jobid}', flush=True)
        except Exception:
            pass  # 展示事件缺席不改变调度器处理结果。
        return result

    executor.handle_job_success = report


def duration(seconds: float) -> str:
    seconds = max(0, int(seconds))
    return f"{seconds // 3600:02d}:{seconds // 60 % 60:02d}:{seconds % 60:02d}"


class LogCursor:
    """增量读取普通行与原生回车进度，限制单次读取和未完成行的内存。"""

    def __init__(self, path: Path):
        self.path = path
        self.offset = 0
        self.identity = None
        self.pending = b""
        self.read_failed = False

    def read(self, final: bool = False) -> list[str]:
        try:
            with self.path.open("rb") as handle:
                st = os.fstat(handle.fileno())
                identity = (st.st_dev, st.st_ino)
                if identity != self.identity or st.st_size < self.offset:
                    self.offset, self.pending = 0, b""
                self.identity = identity
                handle.seek(self.offset)
                chunk = handle.read(256 * 1024)
                self.offset = handle.tell()
                at_end = self.offset == st.st_size
                self.read_failed = False
        except OSError:
            self.read_failed = True
            return []
        fragments = re.split(rb"[\r\n]", self.pending + chunk)
        self.pending = fragments.pop()[-65536:]
        if final and at_end and self.pending:
            fragments.append(self.pending)
            self.pending = b""
        return [ANSI.sub("", value.decode("utf-8", "replace")) for value in fragments if value]


@dataclass
class RunView:
    mode: str
    title: str
    log_path: Path
    project: Path | None = None
    budget: int | None = None
    started: float = field(default_factory=time.monotonic)
    start_text: str = field(default_factory=lambda: datetime.now().astimezone().isoformat(timespec="seconds"))
    status: str = "运行中"
    phase: str = "启动与准备"
    stage: int = 0
    total: int = 0
    phase_started: float = field(default_factory=time.monotonic)
    progress: str = ""
    metadata: str = ""
    result_log: Path | None = None
    detail_log: Path | None = None
    summary_log: Path | None = None
    detail: deque = field(default_factory=lambda: deque(maxlen=100))
    sampling: deque = field(default_factory=lambda: deque(maxlen=100))
    diagnostics: deque = field(default_factory=lambda: deque(maxlen=40))
    warnings: Counter = field(default_factory=Counter)
    phases: dict = field(default_factory=dict)
    routes: dict = field(default_factory=dict)
    jobs: dict = field(default_factory=dict)
    queue_rows: list = field(default_factory=list)
    queue_note: str = ""
    queue_available: bool = False
    resources: dict = field(default_factory=dict)
    failure: str = ""
    cleanup: str = ""
    repaired: bool = False
    worker_success: bool = False
    cells: str = ""
    caches: dict = field(default_factory=dict)
    pending_rule: str = ""
    pending_job: str = ""
    job_context: dict = field(default_factory=dict)
    facts: dict = field(default_factory=dict)
    awaiting_controller_path: bool = False
    awaiting_summary_path: bool = False

    def consume(self, line: str) -> bool:
        """只识别既有文本与明确展示事件；未知原生输出始终保留在完整日志。"""
        line = line.strip()
        if not line:
            return False
        self.diagnostics.append(line[-1200:])
        message = re.sub(r"^\[[^]]+\]\s+(?:INFO|DEBUG|WARNING|ERROR)\s+", "", line)
        significant = False
        for key, token in (('platform','Platform:'),('allocation','检测到当前 Slurm 作业'),
                           ('submission','SLURM job:'),('queue_end','state=COMPLETED exit=0:0'),
                           ('download','下载缺失的 Rust 依赖')):
            if token in message:
                self.facts[key] = message
                significant = True
        submitted = re.match(r'SLURM job:\s*(\d+)', message)
        if submitted:
            self.jobs[submitted[1]] = '构建'
        if self.mode.startswith('doctor-') and re.match(r'(?:创建/更新 |复用未变化的 |在临时目录(?:构建|组装)|构建计划：|下载|编译|control smoke test|其它 Conda 环境创建)', message):
            self.phase = message.split('：',1)[0].split(':',1)[0][:100]
            self.phase_started = time.monotonic()
            self.detail.append(message)
            significant = True
        match = re.search(r"阶段 (\d+)/(\d+)：(.*)|^Step (\d+)/(\d+) START: (.*)", message)
        if match:
            a, b, text = match.groups()[:3] if match[1] else match.groups()[3:]
            self.stage, self.total, self.phase = int(a), int(b), text
            self.phase_started = time.monotonic()
            self.phases[self.stage] = [self.phase, "运行中", "—"]
            significant = True
        done = re.search(r"阶段完成；耗时：(\d+(?:\.\d+)?)s|^Step \d+/\d+ OK: .*?\((\d+)s", message)
        if done and self.stage in self.phases:
            self.phases[self.stage][1:] = ["完成", duration(float(done[1] or done[2]))]
            significant = True
        if message.startswith("提交 QC 表"):
            self.phase, self.phase_started = "最终提交与关闭句柄", time.monotonic()
            significant = True
        if "QC 已完成；" in message:
            self.worker_success = True  # 仍须等待子进程退出（含 finally 关闭句柄）。
        if message.startswith("QC 输入："):
            self.cells = message
        if message.startswith("已复用当前完整 QC"):
            self.caches["CpG"] = "已复用"
        if "单次解码，汇总" in message:
            self.caches["CpG"] = "本次计算"
        if "复用已核验缓存" in message:
            self.caches['CpG'] = '已复用'
        elif '正在从 CpG 表流式提取统计信息' in message:
            self.caches['CpG'] = '本次计算'
        for name in ('Gini', 'TSS'):
            if f'命中 {name}' in message:
                self.caches[name] = '已复用'
            elif f'开始计算 {name}' in message:
                self.caches[name] = '本次计算'
        if '本次跳过 TSS' in message:
            self.caches['TSS'] = '已跳过'
        if "Protocol:" in message:
            self.metadata = message
            significant = True
        if message.startswith("Cores:"):
            try:
                self.budget = int(message.partition(":")[2])
            except ValueError:
                pass
        if "SLURM jobs:" in message or "Mem budget:" in message:
            self.detail.append(message)
            significant = True
        rule = re.match(r'(?:local)?rule ([\w.-]+):', message)
        if rule:
            self.pending_rule, self.pending_job = rule[1], ''
            if self.mode == 'pipeline' and self.phase == '启动与准备':
                self.phase = '工作流调度与执行'
                self.phase_started = time.monotonic()
        local_job = re.match(r'jobid: (\d+)', message)
        if local_job:
            self.pending_job = local_job[1]
            self.job_context[self.pending_job] = [self.pending_rule, '—']
        sample = re.search(r'wildcards:.*?(?:cell|sample|run)=([^, ]+)', message)
        if sample and self.pending_job:
            self.job_context[self.pending_job][1] = sample[1]
        if message.startswith("Controller log:"):
            self.awaiting_controller_path = True
        if self.awaiting_controller_path and self.project is not None:
            # 项目根已通过argv传入；路径可含换行，只有controller文件名是稳定ASCII。
            name = re.search(r'(?:^|/)([A-Za-z0-9_.-]+\.log)$',message)
            if name:
                self.result_log = self.project/'04_logs/controller'/name[1]
                self.awaiting_controller_path = False
        if message.startswith("Detail log :"):
            self.detail_log = Path(message.partition(":")[2].strip())
        if message.startswith("Summary log:"):
            self.awaiting_summary_path = True
        if self.awaiting_summary_path:
            name = re.search(r'(?:^|/)(doctor\.[A-Za-z0-9_.-]+\.summary\.log)$',message)
            if name:
                self.summary_log = self.log_path.parent/name[1]
                self.detail_log = self.log_path.parent/name[1].replace('.summary.log','.detail.log')
                self.awaiting_summary_path = False
        route = re.match(r"\[Alopex route\] (\S+) (running|passed|failed|verified)(?: (.*))?$", message)
        if route:
            name, state, note = route.groups()
            old = self.routes.get(name, {})
            self.routes[name] = {"state": state, "started": old.get("started", time.monotonic()), "note": note or ""}
            significant = True
        if "[REPAIR]" in message:
            self.repaired = True
        if message.startswith("临时沙箱："):
            self.cleanup = message.partition("：")[2].strip()
        native_progress = re.match(r"(\d+) of (\d+) steps \((\d+)%\)", message)
        if native_progress:
            self.progress = f"当前 DAG 已完成 {native_progress[1]}/{native_progress[2]} 个任务"
            significant = True
        job = re.search(rf"Job ({SNAKEMAKE_JOB_ID}) has been submitted with SLURM jobid (\d+(?:_\d+)?)", message)
        if job:
            self.jobs[job[2]] = job[1]
        # Slurm arrays 由官方插件逐元素输出 external_jobid，不能扩成用户全部作业。
        array_job = re.search(rf"external_jobid=(\d+_\d+),\s*snakemake_jobid=({SNAKEMAKE_JOB_ID})(?=\W|$)", message)
        if array_job:
            self.jobs[array_job[1]] = array_job[2]
        finished_job = re.match(rf"(?:Finished jobid:|\[Alopex group finished\]) ({SNAKEMAKE_JOB_ID})(?=\s|$)", message)
        if finished_job:
            task = finished_job[1]
            self.job_context.pop(task, None)
            if self.pending_job == task:
                self.pending_job = ''
            # 采用 controller 的完成事件。squeue 中缺席的任务仍保留给后续观测。
            self.jobs = {job: local for job, local in self.jobs.items() if local != task}
            self.queue_rows = [row for row in self.queue_rows if row['id'] in self.jobs]
        doctor_job = re.match(r"\[Alopex job\] (\d+) (.*)", message)
        if doctor_job:
            self.jobs[doctor_job[1]] = doctor_job[2]
        if "TSS 按 CloneID" in message or message.startswith("CloneID="):
            self.sampling.append(message)
        elif any(word in message for word in ("缓存", "SingleCpG", "开始计算", "跳过 TSS")):
            self.detail.append(message)
        if "WARNING" in line or message.startswith("WARNING:"):
            if message in self.warnings or len(self.warnings) < 100:
                self.warnings[message] += 1
            significant = True
        if re.search(r"\b(?:ERROR|Error|Exception|FAIL)\b|(?:Error|Exception):", line) or "[FAIL]" in line:
            self.failure = message[-1200:]
            significant = True
        return significant

    def finish(self, code: int, interrupted: bool = False) -> None:
        self.status = "已中断" if interrupted else ("READY" if self.mode == "doctor" else "已完成") if code == 0 else "失败"
        if code:
            self.failure = self.failure or f"子进程退出码 {code}；查看完整日志中的错误上下文"
        if self.stage in self.phases and code:
            self.phases[self.stage][1] = self.status

    def lines(self) -> list[str]:
        elapsed = duration(time.monotonic() - self.started)
        finished = self.status != "运行中"
        rows = [f"{self.title}  ·  {self.status}", f"开始  {self.start_text}"]
        if self.project:
            rows.append(f"项目  {self.project.name}")
        if not finished:
            phase = f"{self.stage}/{self.total} · " if self.total else ""
            rows += [f"当前  {phase}{self.phase}", f"耗时  总计 {elapsed} · 本阶段 {duration(time.monotonic() - self.phase_started)}"]
        else:
            rows += [f"总耗时  {elapsed}"]
            if self.resources.get('peak') is not None:
                label={'root_only':'主进程 RSS', 'partial_tree':'可读进程 RSS 合计'}.get(self.resources.get('scope'), '进程 RSS 合计')
                rows.append(f"资源  {label}的采样峰值 {self.resources['peak']/2**30:.2f} GiB")
        if self.metadata:
            rows.append(self.metadata)
        rows.extend(self.facts.values())
        if self.cells:
            rows.append(self.cells)
        if self.caches:
            rows.append('缓存  ' + ' · '.join(f'{name} {state}' for name,state in self.caches.items()))
        if self.progress:
            rows.append("进度  " + self.progress + " · 剩余时间暂无法估计")
        if self.mode.startswith("doctor"):
            for name, item in list(self.routes.items())[:3]:
                label = {"running": "运行中", "passed": "执行完成，待校验", "failed": "失败", "verified": "已通过"}[item['state']]
                rows.append(f"路线  {name} · {label} {item['note']}")
            if finished and self.status == "READY":
                rows.append(f"检查  {self.stage}/{self.total} 阶段通过" if self.total else "检查完成")
                rows.append("修复  " + ("已修复并复验" if self.repaired else "无需修复"))
            if self.cleanup:
                rows.append("清理  " + self.cleanup)
        if not finished:
            r = self.resources
            cpu = f"{r['cpu']:.1f} 核" if r.get('cpu') is not None else "待采集／不可用"
            budget = f" / 预算 {self.budget} 核" if self.budget else ""
            rows.append(f"本机  CPU {cpu}{budget}")
            if r.get('scope') == 'root_only':
                rows.append('采样范围  仅主进程；子树枚举不可用')
            elif r.get('scope') == 'partial_tree':
                rows.append('采样范围  部分进程不可读；CPU/RSS 仅为可读进程合计')
            if r.get('rss') is not None:
                label={'root_only':'主进程 RSS', 'partial_tree':'可读进程 RSS 合计'}.get(r.get('scope'), '进程 RSS 合计')
                rows.append(f"内存  {label} {r['rss']/2**30:.2f} GiB · 采样峰值 {r['peak']/2**30:.2f} GiB")
            else:
                rows.append("内存  进程 RSS 合计不可用")
            count=r.get('children') if r.get('children') is not None else '不可用'
            rows.append(f"进程  主 PID {r.get('pid', '—')} · 子进程 {count} · 活动：{r.get('active', '—')}")
            if self.jobs:
                running = [j for j in self.queue_rows if j['state'] in {'RUNNING', 'COMPLETING'}]
                pending = [j for j in self.queue_rows if j['state'] == 'PENDING']
                counts = f'运行 {len(running)} · 排队 {len(pending)}' if self.queue_available else '运行／排队数量不可用'
                rows.append(f"Slurm  {counts} · {self.queue_note or '等待队列观测'}")
                if self.queue_available:
                    memory = f"{sum(j['memory_mib'] for j in running)/1024:.1f} GiB" if all(j['memory_mib'] is not None for j in running) else '内存不可用'
                    rows.append(f"申请  运行作业 {sum(j['cpus'] for j in running)} 核 / {memory}")
                rows.append("远程实测  CPU/RSS 不可用（不以申请值代替）")
                for item in self.queue_rows[:3]:
                    context = self.job_context.get(self.jobs.get(item['id']), ['阶段未提供','样本未提供'])
                    rows.append(f"Job {item['id']} · {' / '.join(context)} · {item['state']} · {item['runtime']} · {item['node']} · {item['reason']}")
        for warning, count in list(self.warnings.items())[:3]:
            rows.append(f"注意  {warning}" + (f"（重复 {count} 次）" if count > 1 else ""))
        if self.failure and self.status in {"失败", "已中断"}:
            rows += [f"原因  {self.failure}", "下一步  打开完整日志，检查首条错误及对应输入。"]
        if finished and self.status == "已完成" and self.project:
            if self.mode == "qc":
                rows.append(f"输出  {self.project}/06_downstream/DNAme_QC_Information.csv")
            else:
                rows += [f"结果  {self.project}/03_results/", f"报告  {self.project}/03_results/QC_Results/multiqc_report.html"]
        rows.append(f"日志  {self.detail_log or self.result_log or self.log_path}")
        # 控制字符只在展示时转义，完整日志与实际传给OS的路径保持原字节语义。
        return [json.dumps(row,ensure_ascii=False) if any(ord(c)<32 or ord(c)==127 for c in row) else row for row in rows]


class ProcessSampler:
    """采样活跃进程；RSS 求和包含共享页，不用于预算或精确峰值判定。"""

    def __init__(self, pid: int):
        self.pid, self.previous = pid, {}
        self.last_sample = 0.0
        self.root_created = None
        self.peaks = {'tree': 0, 'root_only': 0, 'partial_tree': 0}

    def sample(self, now: float) -> dict:
        try:
            import psutil
            root = psutil.Process(self.pid)
            created = root.create_time()
            if self.root_created is not None and created != self.root_created:
                return dict(pid=self.pid, cpu=None, rss=None, active="原进程已退出")
            self.root_created = created
            try:
                children = root.children(recursive=True)
                tree_available = True
            except (psutil.Error, OSError):
                children, tree_available = [], False
            # 监控目标就是当前展示进程启动的子进程，不能把自有子树一起排除。
            excluded = {os.getpid()}
            found, current, usage, known, rss, active = [], {}, 0.0, 0, 0, []
            incomplete = False
            for proc in [root, *children]:
                if proc.pid in excluded:
                    continue
                try:
                    with proc.oneshot():
                        identity = (proc.pid, proc.create_time())
                        times = proc.cpu_times()
                        cpu_time = times.user + times.system
                        memory = proc.memory_info().rss
                        name = proc.name()
                    current[identity] = cpu_time
                    old = self.previous.get(identity)
                    delta = max(0.0, cpu_time - old) if old is not None else 0.0
                    if old is not None:
                        known += 1
                        usage += delta
                    rss += memory
                    active.append((delta, memory, name))
                    found.append(proc.pid)
                except (psutil.Error, OSError):
                    incomplete = True
                    continue
            interval = now - self.last_sample
            cpu = usage / interval if known and self.last_sample and interval > 0 else None
            self.previous, self.last_sample = current, now
            scope=('partial_tree' if incomplete else 'tree') if tree_available else 'root_only'
            self.peaks[scope] = max(self.peaks[scope], rss)
            return dict(pid=self.pid, cpu=cpu, rss=rss if found else None,
                        peak=self.peaks[scope] if found else None, scope=scope,
                        children=len(children) if found and tree_available else None,
                        active=", ".join(p[2] for p in sorted(active, reverse=True)[:3]))
        except Exception:
            # 可观测性失败不能成为新的科学执行 gate。
            return dict(pid=self.pid, cpu=None, rss=None, active="不可用")


def memory_mib(value: str, cpus: int, nodes: int) -> float:
    match = re.fullmatch(r"(\d+(?:\.\d+)?)([KMGTP]?)([cn]?)", value)
    if not match:
        raise ValueError("unknown Slurm memory unit")
    amount, unit, per = match.groups()
    scale = {'':1, 'K':1/1024, 'M':1, 'G':1024, 'T':1024**2, 'P':1024**3}[unit]
    return float(amount)*scale*(cpus if per == 'c' else nodes)


class SlurmSampler:
    def __init__(self):
        self.next_poll = 0.0

    def sample(self, view: RunView, now: float) -> None:
        if not view.jobs or now < self.next_poll:
            return
        self.next_poll = now + 60
        ids = sorted(view.jobs)
        if any(not re.fullmatch(r"\d+(?:_\d+)?", value) for value in ids):
            return
        try:
            result = subprocess.run(['squeue','--noheader','--array','--jobs',','.join(ids),
                '--format=%i|%T|%M|%N|%r|%C|%m|%D'],capture_output=True,text=True,timeout=5)
            if result.returncode:
                raise OSError(result.stderr.strip())
            rows = []
            for line in result.stdout.splitlines():
                job, state, runtime, node, reason, cpus, memory, nodes = line.split('|')
                if job not in view.jobs:
                    continue
                cpu_count, node_count = int(cpus), int(nodes)
                try:
                    requested_memory = memory_mib(memory,cpu_count,node_count)
                except ValueError:
                    requested_memory = None
                rows.append(dict(id=job,state=state,runtime=runtime,node=node,reason='' if reason=='None' else reason,
                    cpus=cpu_count,memory_mib=requested_memory))
            view.queue_rows = rows
            view.queue_available = True
            view.queue_note = "队列采样 " + datetime.now().astimezone().strftime('%H:%M:%S')
        except (OSError, ValueError, subprocess.SubprocessError):
            self.next_poll = now + 120
            view.queue_note = "队列查询不可用；保留上次观测，不判定任务完成"


class Renderer:
    """一个输出区域；渲染失败切换普通文本，绝不全局接管 stdout/stderr。"""

    def __init__(self, view: RunView, notebook: bool = False, stream: TextIO | None = None):
        self.view, self.stream = view, stream or sys.stdout
        self.live, self.handle = None, None
        self.last_plain, self.last_time = None, 0.0
        try:
            if notebook:
                from IPython import get_ipython
                from IPython.display import display, HTML
                if getattr(get_ipython(), 'kernel', None) is not None:
                    self.handle = display(HTML(self.html()), display_id=True)
            elif self.stream.isatty() and os.environ.get('TERM') != 'dumb':
                from rich.console import Console
                from rich.live import Live
                from rich.text import Text
                self.live = Live(Text(''), console=Console(file=self.stream), auto_refresh=False,
                    screen=False, transient=False, redirect_stdout=False, redirect_stderr=False)
                self.live.start(refresh=False)
        except Exception:
            self.live, self.handle = None, None

    def html(self) -> str:
        view = self.view
        visible = '<br>'.join(html.escape(line) for line in view.lines())
        def details(title, lines, opened=False):
            text = '\n'.join(str(line) for line in lines)
            return f'<details{" open" if opened else ""}><summary>{title}</summary><pre style="white-space:pre-wrap">{html.escape(text)}</pre></details>'
        phases = [' · '.join([str(k), *value]) for k,value in view.phases.items()]
        return ('<div style="border:1px solid #8888;border-radius:8px;padding:12px;line-height:1.65">'
            +visible+details('阶段明细与缓存使用情况', [*phases,*view.detail])
            +details('TSS 抽样明细（完整记录见日志）',view.sampling)
            +details('错误详情' if view.status=='失败' else '日志与诊断（最近40条）',view.diagnostics,view.status=='失败')+'</div>')

    def update(self, significant: bool = False, final: bool = False) -> None:
        now = time.monotonic()
        if not final and now-self.last_time < 1:
            return
        try:
            if self.handle is not None:
                from IPython.display import HTML
                self.handle.update(HTML(self.html()))
            elif self.live is not None:
                from rich.text import Text
                self.live.update(Text('\n'.join(self.view.lines())), refresh=True)
            else:
                key=(self.view.stage,self.view.phase,self.view.status,tuple(self.view.warnings),
                     tuple(self.view.facts.items()),tuple(sorted((r['id'],r['state'],r['reason']) for r in self.view.queue_rows)))
                if not final and key == self.last_plain and now-self.last_time < 300:
                    return
                # 普通日志只写变化和心跳，不逐秒重放面板/路径。
                if self.last_plain is None or final:
                    text='\n'.join(self.view.lines())
                else:
                    text=f'[{datetime.now().astimezone().strftime("%H:%M:%S%z")}] {self.view.status} · {self.view.phase} · 已用 {duration(now-self.view.started)}'
                    if self.view.progress:
                        text+=' · '+self.view.progress
                    if self.view.warnings:
                        text+=' · 注意：'+next(reversed(self.view.warnings))
                    if self.view.facts:
                        text+=' · '+list(self.view.facts.values())[-1]
                print(text,file=self.stream,flush=True)
                self.last_plain=key
            self.last_time=now
        except Exception:
            self._stop_live()
            self.handle=None
            try:
                print(f'显示已退化为普通文本；完整日志：{self.view.log_path}',file=self.stream,flush=True)
            except Exception:
                pass

    def _stop_live(self):
        if self.live is not None:
            try:
                self.live.stop()
            except Exception:
                pass
            self.live=None

    def close(self):
        self._stop_live()


class _Interrupted(BaseException):
    def __init__(self, signum):
        self.signum=signum


def _stop_command(process: subprocess.Popen, signum: int) -> None:
    """只信号本次新建的进程组；等待原入口的清理，超时后回收该组。"""
    try:
        os.killpg(process.pid, signum)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=20)
    except subprocess.TimeoutExpired:
        pass
    # 父命令完成清理或超时后，回收仍在本组运行的后代；wait只等待直接子进程。
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def run_displayed(command: list[str], *, mode: str, title: str, log_dir: Path,
                  project: Path | None = None, budget: int | None = None,
                  env: dict | None = None, notebook: bool = False) -> int:
    """父进程观测独立子进程；输出直接落盘，无未消费 PIPE 或科学文件扫描。"""
    environment=dict(os.environ if env is None else env)
    environment[CHILD_ENV]='1'
    environment['PYTHONUNBUFFERED']='1'
    stamp=datetime.now().astimezone().strftime('%Y%m%dT%H%M%S.%f%z')
    log_path=log_dir/f'{mode}.{stamp}.{os.getpid()}.log'
    log_stream=None
    try:
        log_dir.mkdir(parents=True, exist_ok=True)
        log_stream=log_path.open('xb')
    except OSError as exc:
        try:
            print(f'无法创建展示日志，使用原生输出：{exc}',file=sys.stderr,flush=True)
        except Exception:
            pass
    view=RunView(mode,title,log_path,project,budget) if log_stream is not None else None
    renderer=Renderer(view,notebook=notebook) if view is not None else None
    process=None
    saved={}
    interrupted_by=None

    def interrupt(signum, _frame):
        nonlocal interrupted_by
        if interrupted_by is None:
            interrupted_by=signum
            # 启动期间先记录信号。子进程句柄归属明确后才能转发和回收。
            if process is not None:
                raise _Interrupted(signum)

    try:
        with log_stream if log_stream is not None else nullcontext() as log:
            if threading.current_thread() is threading.main_thread():
                signals=(signal.SIGINT,) if notebook else (signal.SIGINT,signal.SIGTERM,signal.SIGHUP)
                for signum in signals:
                    saved[signum]=signal.signal(signum,interrupt)
            if interrupted_by is not None:
                raise _Interrupted(interrupted_by)
            process=subprocess.Popen(command,env=environment,stdout=log,
                                     stderr=subprocess.STDOUT if log is not None else None,
                                     start_new_session=True)
            if interrupted_by is not None:
                raise _Interrupted(interrupted_by)
            if notebook and saved:
                signal.signal(signal.SIGINT,saved[signal.SIGINT])
                saved.clear()
            # 日志不可用时保留原生流，执行与中断仍走同一进程组所有权边界。
            if log is None:
                code=process.wait()
                return 128-code if code < 0 else code
            sampler=ProcessSampler(process.pid)
            queue=SlurmSampler()
            cursors={log_path:LogCursor(log_path)}
            last_sample=0.0
            while True:
                code=process.poll()
                changed=False
                # 单轮读取有界；结束时仍排空已经落盘的全部日志。
                for cursor in list(cursors.values()):
                    for line in cursor.read(final=code is not None):
                        # Doctor验证调用把stdout重定向到detail；summary文件才包含完整状态。
                        if (mode == 'doctor' and view.summary_log is not None and cursor.path == log_path
                                and not re.match(r'^(?:ERROR:|Traceback|\w+Error:)',line)):
                            continue
                        try:
                            changed=view.consume(line) or changed
                        except Exception:
                            pass  # 展示解析不改变子进程结果；原始文本仍完整落盘。
                    if cursor.read_failed:
                        view.warnings['部分日志暂不可读；保留原日志及子命令退出码'] = 1
                for source in (view.result_log, view.summary_log):
                    if source is not None and source not in cursors:
                        cursors[source]=LogCursor(source)
                now=time.monotonic()
                if now-last_sample >= 2:
                    resources=sampler.sample(now)
                    if resources.get('rss') is None:
                        # 历史峰值仍有效，但失去观测后不能把旧CPU/RSS冒充当前值。
                        resources['peak']=view.resources.get('peak')
                        resources['scope']=view.resources.get('scope')
                    view.resources=resources
                    last_sample=now
                queue.sample(view,now)
                renderer.update(significant=changed)
                if code is not None:
                    try:
                        # 诊断文件失去读取权限不能让已退出的命令永久等待排空。
                        drained=all(c.offset >= c.path.stat().st_size for c in cursors.values()
                                    if not c.read_failed and c.path.exists())
                    except OSError:
                        drained=True
                    if drained:
                        break
                else:
                    time.sleep(0.2)
            code=128-code if code < 0 else code
            view.finish(code)
            if code:
                try:
                    print(f'ERROR: 子命令失败（rc={code}）；日志：{log_path}',file=sys.stderr,flush=True)
                    print('\n'.join(list(view.diagnostics)[-8:]),file=sys.stderr,flush=True)
                except Exception:
                    pass
            try:
                log.write(('\n[Alopex observation] '+json.dumps(dict(exit_code=code,
                    elapsed_seconds=time.monotonic()-view.started,resources=view.resources,
                    resource_scope=('sampled main-process RSS; process-tree enumeration unavailable' if view.resources.get('scope')=='root_only'
                                    else 'sampled readable-process RSS sum only; some processes unavailable'
                                    if view.resources.get('scope')=='partial_tree'
                                    else 'sampled local process-tree RSS sum; shared pages may repeat; not remote or exact peak'
                                    if view.resources.get('scope')=='tree' else 'unavailable')),ensure_ascii=False)+'\n').encode())
            except OSError:
                pass
            renderer.update(final=True)
            return code
    except (KeyboardInterrupt,_Interrupted) as exc:
        signum=exc.signum if isinstance(exc,_Interrupted) else signal.SIGINT
        if process is not None:
            _stop_command(process,signum)
        if view is not None:
            view.finish(128+signum,interrupted=True)
            renderer.update(final=True)
        if notebook:
            raise KeyboardInterrupt from None
        return 128+signum
    finally:
        for signum,handler in saved.items():
            signal.signal(signum,handler)
        if process is not None and process.poll() is None:
            _stop_command(process,signal.SIGTERM)
        if renderer is not None:
            renderer.close()

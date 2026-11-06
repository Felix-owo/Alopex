"""按实际 Python 导入及 Shell 定义追踪阶段的科学源码身份。"""
from __future__ import annotations

import sys

_LOADED_MODULE_CODE = sys._getframe().f_code

import ast
import copy
import hashlib
import importlib.util
import json
import os
import subprocess
import symtable
from collections.abc import Sequence
from functools import lru_cache
from pathlib import Path
from types import CodeType, ModuleType

from .paths import CORE_ROOT


def capture_loaded_source_codes() -> tuple[tuple[str, ModuleType, Path, CodeType], ...]:
    """固定本次已导入的 Alopex 实现；reload 或替换 sys.modules 不能刷新旧的 from-import 引用。"""
    records = []
    for name, module in tuple(sys.modules.items()):
        canonical = getattr(getattr(module, "__spec__", None), "name", None) or name
        if canonical != "alopex" and not canonical.startswith("alopex."):
            continue
        path = Path(module.__file__).resolve()
        if path.name in {"__init__.py", "__main__.py", "plotting.py"}:
            continue
        code = getattr(module, "_LOADED_MODULE_CODE", None)
        if not isinstance(code, CodeType):
            raise RuntimeError(f"Alopex module was loaded before this update; restart the kernel: {name}")
        records.append((name, module, path, code))
    return tuple(records)


def _loaded_constant_structure(value) -> tuple:
    if isinstance(value, CodeType):
        return _loaded_code_structure(value)
    if isinstance(value, (tuple, frozenset)):
        items = tuple(_loaded_constant_structure(item) for item in value)
        return (type(value).__name__, tuple(sorted(items, key=repr)) if isinstance(value, frozenset) else items)
    return type(value).__name__, repr(value)


def _loaded_code_structure(code: CodeType) -> tuple:
    # 行号、源路径和调试行表不定义算法；保留常量的类型/精确字面值及异常表。
    # 只比较 code object，不执行或 reload 磁盘上的模块。
    return (code.co_name, code.co_qualname, code.co_argcount, code.co_posonlyargcount,
            code.co_kwonlyargcount, code.co_nlocals, code.co_stacksize, code.co_flags,
            code.co_code, code.co_exceptiontable, code.co_names, code.co_varnames,
            code.co_freevars, code.co_cellvars,
            tuple(_loaded_constant_structure(value) for value in code.co_consts))



def _loaded_code_contains(module_code: CodeType, active_code: CodeType) -> bool:
    return module_code is active_code or any(
        isinstance(value, CodeType) and _loaded_code_contains(value, active_code)
        for value in module_code.co_consts)


_LOADED_SOURCE_CHECKS: dict[tuple[Path, int], tuple[CodeType, str]] = {}


def _verify_loaded_source(path: Path, source: str | bytes, code: CodeType) -> None:
    data = source.encode("utf-8") if isinstance(source, str) else source
    digest = hashlib.sha256(data).hexdigest()
    key = (path, id(code))
    checked = _LOADED_SOURCE_CHECKS.get(key)
    if checked is not None and checked[0] is code and checked[1] == digest:
        return
    current = compile(source, str(path), "exec", dont_inherit=True, optimize=sys.flags.optimize)
    if _loaded_code_structure(current) != _loaded_code_structure(code):
        raise RuntimeError(f"Alopex source changed after import; restart the kernel: {path}")
    _LOADED_SOURCE_CHECKS[key] = (code, digest)


def assert_loaded_sources_current(records: Sequence[tuple[str, ModuleType, Path, CodeType]]) -> None:
    """拒绝旧内核实现与新磁盘源码混用；成功内容可复用编译比较，逐次读取内容确认未变。"""
    # Notebook 的 from-import 可在 reload 后继续持有旧函数；其 globals 已换新，
    # 因而还须核实际调用栈中的 code object，不能只看 sys.modules 的新记录。
    modules = {id(module.__dict__): code for _, module, _, code in records}
    frame = sys._getframe(1)
    try:
        while frame is not None:
            owner = modules.get(id(frame.f_globals))
            if owner is not None and not _loaded_code_contains(owner, frame.f_code):
                raise RuntimeError("Alopex caller predates a module reload; restart the kernel")
            frame = frame.f_back
    finally:
        del frame
    for name, module, path, code in records:
        if (sys.modules.get(name) is not module
                or getattr(module, "_LOADED_MODULE_CODE", None) is not code):
            raise RuntimeError(f"Alopex was reloaded during this session; restart the kernel: {name}")
        try:
            _verify_loaded_source(path, path.read_bytes(), code)
        except (OSError, SyntaxError) as exc:
            raise RuntimeError(f"Cannot verify loaded Alopex source; restart the kernel: {path}") from exc


def python_ast_bindings(source: str, *, conditional_bindings: bool = False) -> dict[str, ast.AST]:
    """提取去 docstring 的模块绑定；可选保留条件定义的外围条件，不登记推导式局部变量。"""
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
                node.body.pop(0)
    bindings = {}
    def bind(node: ast.AST, owner: ast.AST | None = None) -> None:
        owner = node if owner is None else owner
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            bindings[node.name] = owner
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                for child in ast.walk(target):
                    if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store):
                        bindings[child.id] = owner
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                binding = owner
                if owner is node:
                    # 同一 import 语句内未消费的其他名字不应扩大阶段身份。
                    binding = copy.copy(node)
                    binding.names = [alias]
                bindings[alias.asname or alias.name.split(".")[0]] = binding
        elif conditional_bindings:
            for child in ast.iter_child_nodes(node):
                bind(child, owner)

    for node in tree.body:
        bind(node)
    return bindings


def _global_loads(node: ast.AST) -> set[str]:
    """区分真实全局读取与参数、推导式和嵌套作用域中的局部同名变量。"""
    scope_node = node
    # 去掉唯一 docstring 后的空定义才需要补 Pass，普通绑定保持只读，避免深拷贝。
    if any(isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
           and not child.body for child in ast.walk(node)):
        scope_node = copy.deepcopy(node)
        for child in ast.walk(scope_node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and not child.body:
                child.body.append(ast.Pass())
    table = symtable.symtable(ast.unparse(scope_node), "<scientific binding>", "exec")

    def collect(scope: symtable.SymbolTable) -> set[str]:
        names = {symbol.get_name() for symbol in scope.get_symbols()
                 if symbol.is_global() and symbol.is_referenced()}
        for nested in scope.get_children():
            names.update(collect(nested))
        return names

    return collect(table)


class PythonSourceGraph:
    """只读源码，不 import 或执行被审查模块；每次收集固定一份解析视图。"""

    def __init__(self, source_roots: Sequence[Path] | None = None, *,
                 loaded_sources: Sequence[tuple[str, ModuleType, Path, CodeType]] = ()):
        if source_roots is None:
            source_roots = (CORE_ROOT,)
        self.source_roots = tuple(Path(root).resolve() for root in source_roots)
        for root in self.source_roots:
            if not root.is_dir():
                raise FileNotFoundError(f"Scientific source root is not a directory: {root}")
        self.modules: dict[str, tuple[dict[str, ast.AST], dict[str, tuple[str, str | None]]]] = {}
        self.loaded_codes = {path: code for _, _, path, code in loaded_sources}

    def module_path(self, module: str) -> Path | None:
        if not module or any(not part.isidentifier() for part in module.split(".")):
            raise ValueError(f"Invalid scientific module: {module!r}")
        relative = Path(*module.split("."))
        for root in self.source_roots:
            for candidate in (root / relative.with_suffix(".py"), root / relative / "__init__.py"):
                if candidate.is_file():
                    return candidate
        return None

    def local_namespace(self, module: str) -> bool:
        first = module.split(".")[0]
        return any((root / first).is_dir() or (root / (first + ".py")).is_file()
                   for root in self.source_roots)

    def import_targets(self, module: str, node: ast.AST) -> dict[str, tuple[str, str | None]]:
        if isinstance(node, ast.Import):
            return {alias.asname or alias.name.split(".")[0]: (
                alias.name if alias.asname else alias.name.split(".")[0], None,
            ) for alias in node.names}
        if not isinstance(node, ast.ImportFrom):
            return {}
        path = self.module_path(module)
        package = module if path.name == "__init__.py" else module.rpartition(".")[0]
        if node.level:
            if not package:
                raise ValueError(f"Relative scientific import outside package: {module}")
            target = importlib.util.resolve_name("." * node.level + (node.module or ""), package)
        else:
            target = node.module or ""
        if any(alias.name == "*" for alias in node.names) and self.local_namespace(target):
            raise ValueError(f"Wildcard local scientific import: {module} -> {target}")
        return {alias.asname or alias.name: (target, alias.name) for alias in node.names}

    def load(self, module: str):
        if module in self.modules:
            return self.modules[module]
        path = self.module_path(module)
        if path is None:
            raise ValueError(f"Missing local scientific module: {module}")
        source = path.read_text(encoding="utf-8")
        code = self.loaded_codes.get(path.resolve())
        if code is not None:
            # 校验将要进入 AST 指纹的同一份内容，不能另读已回滚的磁盘版本批准它。
            _verify_loaded_source(path.resolve(), source, code)
        bindings = python_ast_bindings(source, conditional_bindings=True)
        imports: dict[str, tuple[str, str | None]] = {}
        def register(node: ast.AST) -> None:
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                imports.update(self.import_targets(module, node))
            elif not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                for child in ast.iter_child_nodes(node):
                    register(child)

        for node in ast.parse(source).body:
            register(node)
        self.modules[module] = bindings, imports
        return bindings, imports

    def member(self, module: str, parts: Sequence[str]) -> tuple[str, str] | None:
        """解析模块属性链；类/值属性由其完整绑定覆盖，不猜动态模块引用。"""
        if not parts:
            return None
        bindings, imports = self.load(module)
        first, *rest = parts
        if first in imports and rest:
            target, imported_name = imports[first]
            if self.local_namespace(target):
                if imported_name is None:
                    return self.member(target, rest)
                target_bindings, _ = self.load(target)
                if imported_name not in target_bindings:
                    return self.member(target + "." + imported_name, rest)
        if first in bindings:
            return module, first
        child_module = module + "." + first
        if self.module_path(child_module) is not None:
            return self.member(child_module, rest)
        raise ValueError(f"Unknown local scientific dependency: {module}.{'.'.join(parts)}")

    def dependencies(self, roots: Sequence[str], exclude: Sequence[str] = ()) -> dict[str, str]:
        """从本次已读源码视图收集各根的传递依赖；可供多阶段共用解析结果。"""
        selected: dict[str, str] = {}
        excluded = set(exclude)
        pending = []
        for root in roots:
            module, dot, name = root.rpartition(".")
            if not dot:
                raise ValueError(f"Scientific Python root must be module-qualified: {root!r}")
            pending.append((module, name))

        while pending:
            module, name = pending.pop()
            qualified = module + "." + name
            if qualified in selected or qualified in excluded:
                continue
            bindings, imports = self.load(module)
            if name not in bindings:
                raise ValueError(f"Unknown scientific Python root or dependency: {qualified}")
            node = bindings[name]
            selected[qualified] = ast.dump(node, include_attributes=False)

            if name in imports:
                if not isinstance(node, (ast.Import, ast.ImportFrom)):
                    pending.extend((module, dependency) for dependency in
                                   (_global_loads(node) & bindings.keys()) - {name})
                target, imported_name = imports[name]
                if not self.local_namespace(target):
                    continue
                if self.module_path(target) is None:
                    raise ValueError(f"Missing local scientific import: {module} -> {target}")
                if imported_name is not None:
                    dependency = self.member(target, (imported_name,))
                    if dependency is not None:
                        pending.append(dependency)
                continue

            loads = _global_loads(node)
            for dependency in loads & bindings.keys():
                pending.append((module, dependency))

            # 模块属性和函数内惰性导入都须追到真实定义；动态模块分派无法
            # 静态枚举科学消费者，明确失败，不能悄悄留下一个 import 指纹。
            parents = {child: parent for parent in ast.walk(node) for child in ast.iter_child_nodes(parent)}
            all_loads = {child.id for child in ast.walk(node)
                         if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Load)}
            import_uses = [(alias, *imports[alias]) for alias in loads & imports.keys()]
            for child in ast.walk(node):
                if isinstance(child, (ast.Import, ast.ImportFrom)):
                    import_uses.extend((alias, *target) for alias, target in self.import_targets(module, child).items()
                                       if alias in all_loads)
            for alias, target, imported_name in import_uses:
                if not self.local_namespace(target):
                    continue
                target_bindings, _ = self.load(target)
                if imported_name is not None:
                    if imported_name in target_bindings:
                        pending.append((target, imported_name))
                        continue
                    target += "." + imported_name
                self.load(target)
                for child in ast.walk(node):
                    if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Load) and child.id == alias:
                        parent = parents.get(child)
                        if not isinstance(parent, ast.Attribute) or parent.value is not child:
                            raise ValueError(f"Unresolved dynamic local module dependency: {module}.{alias}")
                    if not isinstance(child, ast.Attribute):
                        continue
                    parts = []
                    base = child
                    while isinstance(base, ast.Attribute):
                        parts.append(base.attr)
                        base = base.value
                    if isinstance(base, ast.Name) and base.id == alias:
                        dependency = self.member(target, list(reversed(parts)))
                        if dependency is not None:
                            pending.append(dependency)
        return selected


@lru_cache(maxsize=None)
def scientific_code_identity(
    python: tuple[str, ...] = (), shell: tuple[str, ...] = (),
    *, source_graph: PythonSourceGraph | None = None,
) -> str:
    """按规则实际消费的 Python 根及 Shell 函数计算科学身份，不执行用户启动文件。"""
    graph = source_graph if source_graph is not None else PythonSourceGraph()
    selected = graph.dependencies(python, exclude=("alopex.runtime.effective_cpu_count",))
    shell_code = ""
    if shell:
        identity_env = {key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "ENV"}}
        shell_code = subprocess.check_output(
            ["bash", "--noprofile", "--norc", "-c", 'source "$1"; shift; declare -f "$@"', "_",
             str(CORE_ROOT / "job_runtime.sh"), *shell], text=True, env=identity_env,
        )
    return hashlib.sha256(json.dumps({"python": selected, "shell": shell_code}, sort_keys=True).encode()).hexdigest()

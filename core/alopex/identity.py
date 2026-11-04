"""按实际 Python 导入及 Shell 定义追踪阶段的科学源码身份。"""
from __future__ import annotations

import ast
import copy
import hashlib
import importlib.util
import json
import os
import subprocess
import symtable
from collections.abc import Mapping, Sequence
from functools import lru_cache
from pathlib import Path

from .paths import CORE_ROOT

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

    def __init__(self, source_roots: Sequence[Path] | None = None):
        if source_roots is None:
            source_roots = (CORE_ROOT,)
        self.source_roots = tuple(Path(root).resolve() for root in source_roots)
        for root in self.source_roots:
            if not root.is_dir():
                raise FileNotFoundError(f"Scientific source root is not a directory: {root}")
        self.modules: dict[str, tuple[dict[str, ast.AST], dict[str, tuple[str, str | None]]]] = {}

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


def python_source_dependencies(
    roots: Sequence[str], *, source_roots: Sequence[Path] | None = None,
    exclude: Sequence[str] = (),
) -> dict[str, str]:
    """沿真实本地导入收集限定名 AST；未知本地根/依赖失败，忽略注释与文档。"""
    return PythonSourceGraph(source_roots).dependencies(roots, exclude)


@lru_cache(maxsize=None)
def scientific_code_identity(python: tuple[str, ...] = (), shell: tuple[str, ...] = ()) -> str:
    """按规则实际消费的 Python 根及 Shell 函数计算科学身份，不执行用户启动文件。"""
    selected = python_source_dependencies(python, exclude=("alopex.runtime.effective_cpu_count",))
    shell_code = ""
    if shell:
        identity_env = {key: value for key, value in os.environ.items() if key not in {"BASH_ENV", "ENV"}}
        shell_code = subprocess.check_output(
            ["bash", "--noprofile", "--norc", "-c", 'source "$1"; shift; declare -f "$@"', "_",
             str(CORE_ROOT / "job_runtime.sh"), *shell], text=True, env=identity_env,
        )
    return hashlib.sha256(json.dumps({"python": selected, "shell": shell_code}, sort_keys=True).encode()).hexdigest()

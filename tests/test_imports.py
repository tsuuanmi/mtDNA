"""
Smoke tests verifying that all src/ modules import cleanly and that the
editable-install / absolute-import setup (no sys.path hacks) works
correctly across the codebase.

These tests do NOT require the package to be installed; they rely on the
pytest ``pythonpath = ["."]`` setting in pyproject.toml so that
``from src.*`` works from the project root.
"""

import ast
import importlib
import pkgutil
import re
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _all_src_modules() -> list[str]:
    """Return sorted list of *every* module/package under src/ (recursive)."""
    src_path = Path("src")
    modules: list[str] = []
    # pkgutil.walk_packages recursively discovers all packages and modules
    # (handles any nesting depth, __init__.py packages, and plain .py modules)
    for _, name, _ in pkgutil.walk_packages([str(src_path)], prefix="src."):
        modules.append(name)
    return sorted(modules)


# ---------------------------------------------------------------------------
# Core: every module under src/ must import cleanly
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("module", _all_src_modules())
def test_all_src_modules_import(module: str) -> None:
    """Every Python module and subpackage under src/ (at any depth) must be importable."""
    importlib.import_module(module)  # raises on any error


# ---------------------------------------------------------------------------
# Discovery: every module must expose at least one public callable
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("module", _all_src_modules())
def test_module_has_callable_or_data(module: str) -> None:
    """
    Every module should expose at least one public name (function, class, or
    non-private data).  This catches completely empty stub modules and accidental
    re-exports that pull in nothing.
    """
    mod = importlib.import_module(module)
    names = [name for name in dir(mod) if not name.startswith("_")]
    assert names, f"{module!r} exports no public names — possible empty stub?"


# ---------------------------------------------------------------------------
# Anti-pattern: no sys.path manipulation in any src/ file
# ---------------------------------------------------------------------------


def test_no_sys_path_hacks() -> None:
    """No file under src/ should manipulate sys.path to add the project root."""
    violations: list[str] = []
    for py_file in Path("src").rglob("*.py"):
        content = py_file.read_text()
        for lineno, line in enumerate(content.splitlines(), start=1):
            stripped = line.strip()
            if any(
                pattern in stripped
                for pattern in [
                    "sys.path.insert",
                    "sys.path.append",
                    "sys.path[0:0]",
                    "sys.path = ",
                    "sys.path.extend",
                ]
            ):
                violations.append(f"{py_file}:{lineno}: {line.strip()}")
    assert not violations, "sys.path manipulation found in:\n" + "\n".join(violations)


# ---------------------------------------------------------------------------
# Anti-pattern: no bare 'import sys' without real sys.* usage
# ---------------------------------------------------------------------------


def test_no_standalone_import_sys() -> None:
    """
    Any ``import sys`` in src/ should be accompanied by actual sys.* usage
    (e.g. sys.exit, sys.argv).  A bare import with no sys.* calls is a
    leftover from a removed sys.path hack.
    """
    real_sys = re.compile(
        r"\bsys\.(exit|argv|version|platform|path|modules|"
        r"stdin|stdout|stderr|getrefcount|getsizeof)\b"
    )
    violations: list[str] = []
    for py_file in Path("src").rglob("*.py"):
        content = py_file.read_text()
        if "import sys" not in content:
            continue
        if not real_sys.search(content):
            violations.append(str(py_file))
    assert not violations, "Bare 'import sys' (no sys.* usage) in:\n" + "\n".join(violations)


# ---------------------------------------------------------------------------
# Dead-code guard: no module-level function in src/core may be unreferenced
# ---------------------------------------------------------------------------


def _module_level_functions(root: Path) -> list[tuple[str, int, str]]:
    """Return (file, lineno, name) for every module-level ``def`` under *root*."""
    found: list[tuple[str, int, str]] = []
    for py_file in root.rglob("*.py"):
        try:
            tree = ast.parse(py_file.read_text())
        except SyntaxError:  # pragma: no cover - defensive
            continue
        found.extend((str(py_file), node.lineno, node.name) for node in tree.body if isinstance(node, ast.FunctionDef))
    return found


def _production_name_usage(roots: list[Path]) -> set[str]:
    """Names referenced anywhere under *roots* via ast.Name or ast.Attribute.

    Import aliases (``from x import y``) are deliberately excluded so that a
    pure re-export does not count as real usage and dead functions stay
    detectable.
    """
    used: set[str] = set()
    for root in roots:
        for py_file in root.rglob("*.py"):
            try:
                tree = ast.parse(py_file.read_text())
            except SyntaxError:  # pragma: no cover - defensive
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.Name):
                    used.add(node.id)
                elif isinstance(node, ast.Attribute):
                    used.add(node.attr)
    return used


# Module-level functions in src/core that are intentionally part of the
# documented public API (see docs/core/*.md and docs/variants.md) even though
# no internal production caller currently invokes them. They are exported for
# external boundaries / convenience entry points and must NOT be flagged as
# dead code. Keep this list in sync with the docs.
_DOCUMENTED_PUBLIC_API: frozenset[str] = frozenset(
    {
        "validate_position",  # docs/core/models.md, docs/variants.md — boundary validator contract
        "merge_samples",  # docs/core/mtdna_merger.md — convenience Sample-object merge API
    }
)


def test_no_dead_module_level_functions_in_core() -> None:
    """No *undocumented* module-level function in ``src/core`` may be unreferenced.

    A function is dead if its name never appears as a ``Name``/``Attribute``
    anywhere under ``src/`` or ``scripts/``. Functions referenced only from
    tests (or only in their own docstring/doctest) are legacy slop and must be
    removed. Documented public-API functions with no internal caller are
    exempt (see ``_DOCUMENTED_PUBLIC_API``); they exist for external callers
    and must not be deleted. Methods and class attributes are out of scope.
    """
    production_roots = [Path("src"), Path("scripts")]
    used = _production_name_usage(production_roots)

    dead = [
        f"{rel}:{lineno}:{name}"
        for rel, lineno, name in _module_level_functions(Path("src/core"))
        if name not in used and name not in _DOCUMENTED_PUBLIC_API
    ]
    assert not dead, (
        "Dead module-level functions in src/core (no production references; "
        "only tests/doctests use them). Remove them and their test-only "
        "references:\n  " + "\n  ".join(dead)
    )

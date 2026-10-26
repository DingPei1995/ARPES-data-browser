"""Checks on the shape of the code itself, not on what it computes.

These exist because of a bug that shipped: the packages are called
``loader``, ``tools`` and ``ui``, and ``loader`` is also the obvious name
for a variable holding one, so

    for loader in loader.registry.loaders():

made ``loader`` local to the function and the right-hand side raised
``UnboundLocalError`` -- at runtime, in the one branch that opened the Load
window. Nothing catches that at import time, and no test of *behaviour*
catches it either unless the test happens to run that branch.

So: walk the ASTs and check the property directly.
"""
import ast
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGES = ("loader", "tools", "ui", "devtools")


def _python_files():
    for folder, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d != "__pycache__"]
        for name in sorted(files):
            if name.endswith(".py"):
                yield os.path.join(folder, name)


def _parse(path):
    with open(path, encoding="utf-8") as handle:
        return ast.parse(handle.read(), filename=path)


def _imported_names(tree):
    """The names a module's imports bind. ``import a.b`` binds ``a``."""
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                names.add((alias.asname or alias.name).split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                names.add(alias.asname or alias.name)
    return names


def _bound_names(scope):
    bound = set()
    for child in ast.walk(scope):
        if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store):
            bound.add(child.id)
        elif isinstance(child, ast.arg):
            bound.add(child.arg)
        elif isinstance(child, ast.ExceptHandler) and child.name:
            bound.add(child.name)
        elif isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if child is not scope:
                bound.add(child.name)
    return bound


def _used_as_module(scope, name):
    for child in ast.walk(scope):
        if (isinstance(child, ast.Attribute)
                and isinstance(child.value, ast.Name)
                and child.value.id == name):
            return True
    return False


def test_no_local_variable_shadows_an_imported_module():
    """The regression this file exists for."""
    problems = []
    for path in _python_files():
        tree = _parse(path)
        imports = _imported_names(tree)
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            declared = {name for child in ast.walk(node)
                        if isinstance(child, (ast.Global, ast.Nonlocal))
                        for name in child.names}
            for name in sorted((_bound_names(node) & imports) - declared):
                if _used_as_module(node, name):
                    problems.append(
                        f"{os.path.relpath(path, ROOT)}:{node.lineno} "
                        f"{node.name}() binds '{name}', which is also an "
                        f"imported module used as one in the same scope")
    assert not problems, "\n".join(problems)


def test_packages_are_never_referenced_by_their_bare_name():
    """``from tools import kspace`` rather than ``import tools.kspace``.

    The second binds ``tools``, so every use reads ``tools.kspace.X`` and
    any variable called ``tools`` shadows it. Binding the module name
    instead makes a collision a plain NameError at import time rather than
    an UnboundLocalError down one branch.

    The launcher is exempt: its ``ui`` is the built main window, a global it
    uses on nearly every line, and it imports the ``ui`` package's modules
    under their own names precisely so the two cannot collide.
    """
    offenders = []
    for path in _python_files():
        tree = _parse(path)
        for node in ast.walk(tree):
            if not isinstance(node, ast.Import):
                continue
            for alias in node.names:
                head = alias.name.split(".")[0]
                if head in PACKAGES and "." in alias.name and not alias.asname:
                    offenders.append(
                        f"{os.path.relpath(path, ROOT)}:{node.lineno} "
                        f"import {alias.name}")
    assert not offenders, "\n".join(offenders)


def test_every_module_imports_cleanly():
    """Nothing in the tree fails on import -- which a typo in a rewritten
    import statement would, and which no other test would notice for a
    module that only gets imported down one menu path."""
    import importlib
    import pkgutil
    import sys

    if ROOT not in sys.path:
        sys.path.insert(0, ROOT)
    failures = []
    for package in PACKAGES:
        for info in pkgutil.iter_modules([os.path.join(ROOT, package)]):
            name = f"{package}.{info.name}"
            if package == "devtools":
                continue          # needs spglib / the lab's .mat file
            try:
                importlib.import_module(name)
            except Exception as exc:                        # noqa: BLE001
                failures.append(f"{name}: {type(exc).__name__}: {exc}")
    assert not failures, "\n".join(failures)


def test_the_qt_free_packages_really_are():
    """``loader`` and ``tools`` are importable on a machine with no display,
    which is what makes them usable from a script -- and is the reason the
    tree is split this way at all. A stray ``from PyQt5 ...`` in either
    would quietly undo that."""
    offenders = []
    for path in _python_files():
        relative = os.path.relpath(path, ROOT)
        if not (relative.startswith("loader" + os.sep)
                or relative.startswith("tools" + os.sep)):
            continue
        tree = _parse(path)
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module]
            for name in names:
                if name.split(".")[0] in ("PyQt5", "pyqtgraph"):
                    offenders.append(f"{relative}:{node.lineno} imports {name}")
    assert not offenders, "\n".join(offenders)


@pytest.mark.parametrize("package", PACKAGES + ("test",))
def test_every_package_says_what_it_is_for(package):
    init = os.path.join(ROOT, package, "__init__.py")
    assert os.path.isfile(init), f"{package} is not a package"
    with open(init, encoding="utf-8") as handle:
        assert ast.get_docstring(ast.parse(handle.read())), \
            f"{package}/__init__.py has no docstring"

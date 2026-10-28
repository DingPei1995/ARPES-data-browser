"""The macOS launcher, ``ARPES viewer.app``: a double-click starts the
program. Checked here without a Mac -- that the bundle is one macOS will
open, that its script is executable (in git too, since a clone or a ZIP
takes the mode from there) and parses, and that its version is the
program's."""
import os
import plistlib
import subprocess

import pytest

import version

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
APP = os.path.join(ROOT, "ARPES viewer.app", "Contents")


def _info():
    with open(os.path.join(APP, "Info.plist"), "rb") as handle:
        return plistlib.load(handle)


def test_the_bundle_names_its_executable_and_icon():
    info = _info()
    assert info["CFBundlePackageType"] == "APPL"
    script = os.path.join(APP, "MacOS", info["CFBundleExecutable"])
    assert os.path.isfile(script)
    assert os.path.isfile(os.path.join(APP, "Resources",
                                       info["CFBundleIconFile"] + ".icns"))


def test_the_bundle_version_is_the_programs():
    info = _info()
    assert info["CFBundleShortVersionString"] == version.VERSION
    assert info["CFBundleVersion"] == version.VERSION


def test_the_script_is_executable_and_parses():
    script = os.path.join(APP, "MacOS", _info()["CFBundleExecutable"])
    with open(script, encoding="utf-8") as handle:
        assert handle.readline().startswith("#!/bin/bash")
    if os.name == "posix":
        assert os.access(script, os.X_OK)
    subprocess.run(["bash", "-n", script], check=True)


def test_the_script_is_executable_in_git():
    try:
        out = subprocess.run(
            ["git", "ls-files", "-s", "ARPES viewer.app/Contents/MacOS"],
            cwd=ROOT, capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("not a git checkout")
    if not out.strip():
        pytest.skip("not committed yet")
    assert out.split()[0] == "100755", out

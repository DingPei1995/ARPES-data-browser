"""Number boxes use a decimal point whatever the system locale is.

Run in a fresh interpreter under a French system locale: Qt reads the
locale once, when it starts, so this test cannot share the suite's process.
"""
import os
import subprocess
import sys

import pytest

pytest.importorskip("PyQt5")

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))

SCRIPT = """
from PyQt5.QtCore import QLocale
from PyQt5.QtWidgets import QApplication, QDoubleSpinBox
assert QLocale.system().decimalPoint() == ",", QLocale.system().name()
import ARPES_viewer  # noqa: F401
app = QApplication([])
box = QDoubleSpinBox()
box.setDecimals(3)
box.setRange(-10000, 10000)
box.setValue(1234.5)
print(box.text(), box.valueFromText("0.05"))
"""


def test_spin_boxes_use_a_decimal_point_under_a_french_locale():
    env = dict(os.environ, LC_ALL="fr_FR.UTF-8", LANG="fr_FR.UTF-8",
               QT_QPA_PLATFORM="offscreen")
    result = subprocess.run([sys.executable, "-c", SCRIPT], cwd=ROOT, env=env,
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert result.stdout.split() == ["1234.500", "0.05"]

"""Tum EDA adimlarini sirayla calistirir (veriler indirildikten sonra).

Kullanim: python scripts/run_eda.py
Onkosul: data/golub/ ve data/armstrong/MLL.zip yerinde; 09 icin internet (GEO platform tablolari).
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEPS = [
    "01_validate_data.py",
    "02_eda_quality.py",
    "06_armstrong_inspect.py",
    "08_armstrong_validate.py",
    "09_common_genes.py",
    "10_eda_projection.py",
    "11_batch_checks.py",
    "12_joint_eda.py",
    "13_missing_and_clustering.py",
]
for step in STEPS:
    print(f"\n>>> {step}")
    subprocess.run([sys.executable, str(ROOT / "scripts" / step)], check=True, cwd=ROOT)
print("\nTum adimlar tamamlandi. Sekiller ve tablolar: figures/")

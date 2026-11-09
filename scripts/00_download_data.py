"""Golub verisini Kaggle'dan indirir (ham veri depoya konmaz).
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config

KAGGLE_SLUG = "crawford/gene-expression"  

config.GOLUB_DIR.mkdir(parents=True, exist_ok=True)
subprocess.run(
    ["kaggle", "datasets", "download", "-d", KAGGLE_SLUG,
     "-p", str(config.GOLUB_DIR), "--unzip"],
    check=True,
)
print("İndirildi ->", config.GOLUB_DIR)

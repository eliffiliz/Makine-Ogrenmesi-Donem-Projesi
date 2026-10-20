"""Golub verisini Kaggle'dan indirir (ham veri depoya konmaz).

Önkoşul: Kaggle hesabı ve ~/.kaggle/kaggle.json (API token).
Çalıştır:  python scripts/00_download_data.py
Kaggle veri seti adını (slug) sayfadan doğrula: "Gene expression dataset (Golub)".
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config

KAGGLE_SLUG = "crawford/gene-expression"  # DOĞRULA: Kaggle sayfasındaki adla eşleşmeli

config.GOLUB_DIR.mkdir(parents=True, exist_ok=True)
subprocess.run(
    ["kaggle", "datasets", "download", "-d", KAGGLE_SLUG,
     "-p", str(config.GOLUB_DIR), "--unzip"],
    check=True,
)
print("İndirildi ->", config.GOLUB_DIR)

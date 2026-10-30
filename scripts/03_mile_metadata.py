"""MILE (GSE13159) ornek bilgilerini (metadata) indirir ve ozetler.

Ifade matrisini YUKLEMEZ; yalnizca '!Sample_*' baslik satirlarini okur.

Kullanim:
  python scripts/03_mile_metadata.py                 # otomatik indirmeyi dener
  python scripts/03_mile_metadata.py DOSYA.txt.gz    # elle indirdigin dosyayi kullan
"""
import gzip
import shutil
import sys
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config

MILE_DIR = config.DATA_DIR / "mile"
URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13159/matrix/"
       "GSE13159_series_matrix.txt.gz")


def get_file(arg):
    if arg:
        return Path(arg)
    MILE_DIR.mkdir(parents=True, exist_ok=True)
    dest = MILE_DIR / "GSE13159_series_matrix.txt.gz"
    if not dest.exists():
        print("Indiriliyor:", URL)
        with urllib.request.urlopen(URL) as r, open(dest, "wb") as f:
            shutil.copyfileobj(r, f)
    return dest


def clean(s):
    return s.strip().strip('"')


def parse_header(path):
    """'!series_matrix_table_begin' satirina kadar olan '!Sample_' satirlarini okur."""
    fields = {}
    sample_ids = None
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                fields.setdefault(parts[0], []).append([clean(p) for p in parts[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                header = next(f).rstrip("\n").split("\t")
                sample_ids = [clean(h) for h in header[1:]]
                break
    return fields, sample_ids


def build_table(fields, sample_ids):
    n = len(sample_ids)
    df = pd.DataFrame({"geo_accession": sample_ids})
    for key in ("!Sample_title", "!Sample_source_name_ch1", "!Sample_platform_id"):
        if key in fields:
            df[key.replace("!Sample_", "")] = fields[key][0][:n]
    for i, row in enumerate(fields.get("!Sample_characteristics_ch1", [])):
        vals = row[:n]
        keys = {v.split(":", 1)[0].strip() for v in vals if ":" in v}
        if len(keys) == 1:
            name = keys.pop()
            df[f"char_{name}"] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
        else:
            df[f"char_line{i}"] = vals
    return df


def main():
    path = get_file(sys.argv[1] if len(sys.argv) > 1 else None)
    fields, sample_ids = parse_header(path)
    df = build_table(fields, sample_ids)

    MILE_DIR.mkdir(parents=True, exist_ok=True)
    out = MILE_DIR / "mile_metadata.csv"
    df.to_csv(out, index=False)

    print(f"Ornek sayisi: {len(df)}")
    print("Sutunlar:", list(df.columns))
    print(f"Kaydedildi -> {out}\n")
    for col in df.columns:
        if col == "geo_accession":
            continue
        vc = df[col].value_counts()
        print(f"--- {col} ({len(vc)} farkli deger) ---")
        print(vc.head(40).to_string())
        print()
    print("Ornek satirlar (ilk 5):")
    print(df.head(5).to_string())


if __name__ == "__main__":
    main()

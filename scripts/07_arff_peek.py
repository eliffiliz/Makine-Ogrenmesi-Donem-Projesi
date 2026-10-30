"""Armstrong MLL.arff dosyasinin yapisini gosterir: ozellik adlari, sinif etiketi, boyut, deger araligi.

Girdi: data/armstrong/raw/MLL.arff   (06_armstrong_inspect.py zip'i acar)
Bu betik yalnizca dosyayi OKUR ve ozetler; veriyi degistirmez.
"""
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config

p = config.DATA_DIR / "armstrong" / "raw" / "MLL.arff"
rel, attrs, rows, in_data = None, [], [], False
with open(p, encoding="utf-8", errors="replace") as f:
    for line in f:
        s = line.strip()
        if not s or s.startswith("%"):
            continue
        low = s.lower()
        if in_data:
            rows.append(s)
        elif low.startswith("@relation"):
            rel = s
        elif low.startswith("@attribute"):
            attrs.append(s)
        elif low.startswith("@data"):
            in_data = True

print("Relation:", rel)
print(f"@attribute sayisi: {len(attrs)}")
print("Ilk 8 ozellik tanimi:")
for a in attrs[:8]:
    print("  ", a[:120])
print("Son 3 ozellik tanimi (sinif etiketi burada olmali):")
for a in attrs[-3:]:
    print("  ", a[:200])

print(f"\nVeri satiri sayisi: {len(rows)}")
if rows and rows[0].startswith("{"):
    print("UYARI: seyrek (sparse) ARFF bicimi, ayri islem gerekir.")
    sys.exit()

split = [r.split(",") for r in rows]
widths = Counter(len(r) for r in split)
print("Satir basina alan sayisi (adet):", dict(widths))
labels = Counter(r[-1].strip() for r in split)
print("Sinif etiketleri (son alan):", dict(labels))

vals = np.array([[float(x) for x in r[:-1]] for r in split], dtype=float)
print(f"\nDeger matrisi: {vals.shape[0]} ornek x {vals.shape[1]} ozellik")
print(f"NaN: {int(np.isnan(vals).sum())}, min = {np.nanmin(vals):.3f}, max = {np.nanmax(vals):.3f}, "
      f"medyan = {np.nanmedian(vals):.3f}, negatif orani = {(vals < 0).mean():.3f}")
print("Ilk satirin ilk 5 degeri:", np.round(vals[0, :5], 3).tolist())
print("Ilk 5 satirin etiketi:", [r[-1].strip() for r in split[:5]])
print("Son 5 satirin etiketi:", [r[-1].strip() for r in split[-5:]])

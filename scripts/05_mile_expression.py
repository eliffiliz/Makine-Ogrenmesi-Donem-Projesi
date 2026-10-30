"""MILE (GSE13159): ALL/AML ornekleri icin ifade matrisini okur, kontrol eder ve .npy olarak onbellege alir.

Girdi : data/mile/GSE13159_series_matrix.txt.gz  ve  data/mile/mile_labels_all_aml.csv (04 uretir)
Cikti : data/mile/mile_expr_all_aml.npy  (ornek x prob, float32)
        data/mile/mile_probes.txt, data/mile/mile_samples.csv
        figures/mile_value_distribution.(png|pdf)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config

MILE = config.DATA_DIR / "mile"
MATRIX = MILE / "GSE13159_series_matrix.txt.gz"

labels = pd.read_csv(MILE / "mile_labels_all_aml.csv")
wanted = set(labels["geo_accession"])
print(f"Okunacak ornek: {len(wanted)} (birkac dakika surebilir, bekle)")

# '!' ile baslayan satirlar (GEO basligi ve tablo sonu) yorum sayilip atlanir
expr = pd.read_csv(MATRIX, sep="\t", comment="!", index_col=0,
                   usecols=lambda c: c == "ID_REF" or c in wanted)
expr = expr.astype(np.float32)

labels = labels.set_index("geo_accession").reindex(expr.columns)
assert labels["label"].notna().all(), "Bazi ornekler etiket tablosunda bulunamadi"

X = expr.T.to_numpy()  # ornek x prob
print(f"Matris: {X.shape[0]} ornek x {X.shape[1]} prob")

v = X.ravel()
print(f"NaN sayisi: {int(np.isnan(v).sum())}")
print(f"min = {np.nanmin(v):.3f}, max = {np.nanmax(v):.3f}")
pcts = [0.1, 1, 5, 25, 50, 75, 95, 99, 99.9]
print("Yuzdelikler:", {p: round(float(np.nanpercentile(v, p)), 3) for p in pcts})
print(f"Negatif oran: {(v < 0).mean():.4f}")
if np.nanmax(v) < 30:
    print("Tahmin: degerler zaten log2 olcekli gorunuyor (maks < 30). Tekrar log alma!")
else:
    print("Tahmin: degerler ham/lineer olcekli gorunuyor (maks >= 30).")

med = np.nanmedian(X, axis=1)
summary = pd.DataFrame({"median": med, "lineage": labels["lineage"].values,
                        "sample_type": labels["sample_type"].values})
print("\nOrnek basina medyan ifade, gruplara gore:")
print(summary.groupby(["lineage", "sample_type"])["median"].agg(["count", "mean", "std"]).round(3).to_string())

np.save(MILE / "mile_expr_all_aml.npy", X)
(MILE / "mile_probes.txt").write_text("\n".join(expr.index.astype(str)), encoding="utf-8")
labels.reset_index().to_csv(MILE / "mile_samples.csv", index=False)
print("\nKaydedildi ->", MILE)

# Sekil: tum degerlerin dagilimi ve ornek medyanlari
config.FIG_DIR.mkdir(parents=True, exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(8, 3))
rng = np.random.default_rng(config.SEED)
axes[0].hist(rng.choice(v, size=min(2_000_000, v.size), replace=False), bins=100, density=True)
axes[0].set_xlabel("Ifade degeri (GEO matrisi)")
axes[0].set_ylabel("Yogunluk")
groups = [g["median"].values for _, g in summary.groupby(["lineage", "sample_type"])]
names = [f"{a}\n{b}" for (a, b), _ in summary.groupby(["lineage", "sample_type"])]
axes[1].boxplot(groups, showmeans=True)
axes[1].set_xticklabels(names, fontsize=7)
axes[1].set_ylabel("Ornek basina medyan ifade")
fig.tight_layout()
fig.savefig(config.FIG_DIR / "mile_value_distribution.png", dpi=300)
fig.savefig(config.FIG_DIR / "mile_value_distribution.pdf")
print("Sekil ->", config.FIG_DIR)

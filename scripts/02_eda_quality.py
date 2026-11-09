"""EDA adim 1: Golub veri kalitesi ve ifade dagilimlari (log oncesi/sonrasi).

Ciktilar (figures/):
  eda_class_distribution.(png|pdf)
  eda_value_distribution.(png|pdf)
  eda_sample_medians.(png|pdf)
  eda_quality_summary.csv

NOT: Alt/ust esik (taban 100, tavan 16000; src/eda.py) literaturdeki sabit sabitlerdir (Dudoit ve ark. 2002),
veriden ogrenilmez, bu yuzden sizinti olusturmaz.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config, eda
from src.data_loading import load_golub

FLOOR, CEIL = eda.FLOOR, eda.CEIL
plt.rcParams.update({"font.size": 10, "figure.dpi": 100})
config.FIG_DIR.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / f"{name}.png", dpi=300)
    fig.savefig(config.FIG_DIR / f"{name}.pdf")
    plt.close(fig)


X_train, y_train, X_test, y_test = load_golub()
sets = {"Training (n=38)": (X_train, y_train), "Independent test (n=34)": (X_test, y_test)}
colors = {"Training (n=38)": "tab:blue", "Independent test (n=34)": "tab:orange"}

# --- 1) Class distribution ---
fig, ax = plt.subplots(figsize=(4.8, 3))
names = list(sets)
all_c = [int((sets[n][1] == 0).sum()) for n in names]
aml_c = [int((sets[n][1] == 1).sum()) for n in names]
x = np.arange(len(names))
ax.bar(x - 0.2, all_c, 0.4, label="ALL")
ax.bar(x + 0.2, aml_c, 0.4, label="AML")
for i, (a, b) in enumerate(zip(all_c, aml_c)):
    ax.text(i - 0.2, a + 0.3, str(a), ha="center")
    ax.text(i + 0.2, b + 0.3, str(b), ha="center")
ax.set_xticks(x)
ax.set_xticklabels(names)
ax.set_ylabel("Number of samples")
ax.legend()
save(fig, "eda_class_distribution")

# --- 2) Raw and log2 distributions ---
fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.2))
rows = []
for name, (X, y) in sets.items():
    v = X.to_numpy(dtype=float).ravel()
    lo, hi = np.percentile(v, [1, 99])
    axes[0].hist(v[(v >= lo) & (v <= hi)], bins=100, density=True, alpha=0.5, color=colors[name], label=name)
    axes[1].hist(np.log2(np.clip(v, FLOOR, CEIL)), bins=100, density=True, alpha=0.5, color=colors[name], label=name)
    rows.append({
        "set": name, "n_values": v.size,
        "negative_ratio": float((v < 0).mean()),
        f"below_{FLOOR:g}_ratio": float((v < FLOOR).mean()),
        f"above_{CEIL:g}_ratio": float((v > CEIL).mean()),
        "min": v.min(), "median": float(np.median(v)), "max": v.max(),
    })
axes[0].set_title("Raw values (1st-99th percentile)")
axes[0].set_xlabel("Expression value")
axes[0].set_ylabel("Density")
axes[1].set_title(f"log2(clip[{FLOOR:g}, {CEIL:g}])")
axes[1].set_xlabel("log2 expression")
axes[1].legend()
save(fig, "eda_value_distribution")

summary = pd.DataFrame(rows)
summary.to_csv(config.FIG_DIR / "eda_quality_summary.csv", index=False)
print(summary.to_string(index=False))

# --- 3) Per-sample median (batch hint) ---
fig, ax = plt.subplots(figsize=(4.8, 3.2))
data, labels = [], []
for name, (X, y) in sets.items():
    data.append(np.median(np.log2(np.clip(X.to_numpy(dtype=float), FLOOR, CEIL)), axis=1))
    labels.append(name)
ax.boxplot(data, showmeans=True)
ax.set_xticklabels(labels)
for i, d in enumerate(data, start=1):
    ax.scatter(np.random.default_rng(config.SEED).normal(i, 0.04, len(d)), d, s=10, alpha=0.6, color="k")
ax.set_ylabel("Per-sample median log2 expression")
save(fig, "eda_sample_medians")
print("Figures ->", config.FIG_DIR)

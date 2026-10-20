"""EDA adım 1: veri kalitesi ve ifade dağılımları (log öncesi/sonrası).

Üretilen çıktılar (figures/):
  eda_class_distribution.(png|pdf)
  eda_value_distribution.(png|pdf)
  eda_sample_medians.(png|pdf)
  eda_quality_summary.csv

NOT: Alt/üst eşik (20, 16000) sabit sabitlerdir, veriden öğrenilmez; bu yüzden
sızıntı oluşturmaz. Bu betik yalnızca keşif amaçlıdır, model eğitmez.
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
from src.data_loading import load_golub

FLOOR, CEIL = 20.0, 16000.0  # rehberdeki örnek eşik aralığı
plt.rcParams.update({"font.size": 10, "figure.dpi": 100})
config.FIG_DIR.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / f"{name}.png", dpi=300)
    fig.savefig(config.FIG_DIR / f"{name}.pdf")
    plt.close(fig)


X_train, y_train, X_test, y_test = load_golub()
sets = {"Eğitim (n=38)": (X_train, y_train), "Bağımsız test (n=34)": (X_test, y_test)}
colors = {"Eğitim (n=38)": "tab:blue", "Bağımsız test (n=34)": "tab:orange"}

# --- 1) Sınıf dağılımı ---------------------------------------------------
fig, ax = plt.subplots(figsize=(4.5, 3))
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
ax.set_ylabel("Örnek sayısı")
ax.legend()
save(fig, "eda_class_distribution")

# --- 2) Ham ve log2 dağılımlar -------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(8, 3))
rows = []
for name, (X, y) in sets.items():
    v = X.to_numpy(dtype=float).ravel()
    lo, hi = np.percentile(v, [1, 99])
    axes[0].hist(v[(v >= lo) & (v <= hi)], bins=100, density=True,
                 alpha=0.5, color=colors[name], label=name)
    v_log = np.log2(np.clip(v, FLOOR, CEIL))
    axes[1].hist(v_log, bins=100, density=True, alpha=0.5,
                 color=colors[name], label=name)
    rows.append({
        "küme": name,
        "n_değer": v.size,
        "negatif_oran": float((v < 0).mean()),
        f"<{FLOOR:g}_oran": float((v < FLOOR).mean()),
        f">{CEIL:g}_oran": float((v > CEIL).mean()),
        "min": v.min(), "medyan": float(np.median(v)), "max": v.max(),
    })
axes[0].set_title("Ham değerler (1.-99. yüzdelik)")
axes[0].set_xlabel("İfade değeri")
axes[0].set_ylabel("Yoğunluk")
axes[1].set_title(f"log2(kırp[{FLOOR:g}, {CEIL:g}])")
axes[1].set_xlabel("log2 ifade")
axes[1].legend()
save(fig, "eda_value_distribution")

summary = pd.DataFrame(rows)
summary.to_csv(config.FIG_DIR / "eda_quality_summary.csv", index=False)
print(summary.to_string(index=False))

# --- 3) Örnek başına medyan (batch ipucu) ---------------------------------
fig, ax = plt.subplots(figsize=(4.5, 3))
data, labels = [], []
for name, (X, y) in sets.items():
    med = np.median(np.log2(np.clip(X.to_numpy(dtype=float), FLOOR, CEIL)), axis=1)
    data.append(med)
    labels.append(name)
ax.boxplot(data, showmeans=True)
ax.set_xticklabels(labels)
for i, d in enumerate(data, start=1):
    ax.scatter(np.random.default_rng(config.SEED).normal(i, 0.04, len(d)), d,
               s=10, alpha=0.6, color="k")
ax.set_ylabel("Örnek başına medyan log2 ifade")
save(fig, "eda_sample_medians")

print("Şekiller ->", config.FIG_DIR)

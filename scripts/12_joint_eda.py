"""EDA adim 4: Armstrong veri kalitesi + Golub-Armstrong birlesik analiz (ortak genler).

A) Armstrong: sinif dagilimi, ham ve log2 dagilim, ornek medyanlari, ozet tablo (A1: veri kalitesi)
B) Birlesik PCA (ortak gen sembolleri): ham birlestirme ile veri seti ici z-skor karsilastirmasi.
   Veri setleri arasi dagilim/batch farkini gosterir (A1: iki veri seti arasi fark).
   Her sembol icin en yuksek varyansli prob secilir (etiket kullanmaz); Armstrong'da MLL -> ALL.
Onkosul: 09_common_genes.py calistirilmis olmali (platform tablolari data/annotation/ icinde).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config, eda
from src.armstrong_loading import binary_labels, load_armstrong
from src.data_loading import load_golub

plt.rcParams.update({"font.size": 9})
config.FIG_DIR.mkdir(parents=True, exist_ok=True)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(config.FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------- A) Armstrong veri kalitesi ----------------
XA, y3 = load_armstrong()
v = XA.to_numpy().ravel()
fig, ax = plt.subplots(1, 4, figsize=(14, 3))
counts = y3.value_counts().reindex(["ALL", "MLL", "AML"])
ax[0].bar(counts.index, counts.values, color=["tab:blue", "tab:purple", "tab:red"])
for i, c in enumerate(counts.values):
    ax[0].text(i, c + 0.4, str(c), ha="center")
ax[0].set_ylabel("Number of samples")
ax[0].set_title("Class distribution", fontsize=9)
lo, hi = np.percentile(v, [1, 99])
ax[1].hist(v[(v >= lo) & (v <= hi)], bins=100, density=True)
ax[1].set_title("Raw values (1st-99th percentile)", fontsize=9)
ax[1].set_xlabel("Expression value")
ax[1].set_ylabel("Density")
ax[2].hist(np.log2(np.clip(v, eda.FLOOR, eda.CEIL)), bins=100, density=True)
ax[2].set_title(f"log2(clip[{eda.FLOOR:g}, {eda.CEIL:g}])", fontsize=9)
ax[2].set_xlabel("log2 expression")
med = np.median(eda.to_log2(XA.to_numpy()), axis=1)
groups = [med[y3.to_numpy() == c] for c in ["ALL", "MLL", "AML"]]
ax[3].boxplot(groups, showmeans=True)
ax[3].set_xticklabels(["ALL", "MLL", "AML"])
ax[3].set_ylabel("Per-sample median log2 expression")
ax[3].set_title("Sample medians by class", fontsize=9)
save(fig, "eda_armstrong_quality")
summ = pd.DataFrame([{
    "dataset": "Armstrong", "n_values": v.size, "negative_ratio": float((v < 0).mean()),
    f"below_{eda.FLOOR:g}_ratio": float((v < eda.FLOOR).mean()),
    f"above_{eda.CEIL:g}_ratio": float((v > eda.CEIL).mean()),
    "min": v.min(), "median": float(np.median(v)), "max": v.max()}])
summ.to_csv(config.FIG_DIR / "eda_armstrong_quality_summary.csv", index=False)
print(summ.round(4).to_string(index=False))

# ---------------- B) Birlesik analiz ----------------
X_tr, y_tr, X_te, y_te = load_golub()
raw_ids = pd.read_csv(config.TRAIN_FILE, usecols=["Gene Accession Number"])["Gene Accession Number"].astype(str).to_numpy()
Xg_all = pd.concat([X_tr, X_te]).to_numpy()
yg = pd.concat([y_tr, y_te]).to_numpy()
ctrl = np.array([i.startswith("AFFX-") for i in raw_ids])
Xg = Xg_all[:, ~ctrl]
idg = raw_ids[~ctrl]
dmg = eda.dudoit_mask(Xg)
Lg, idg = eda.to_log2(Xg)[:, dmg], idg[dmg]
symg = eda.probe_symbols("GPL80", idg)

XAv = XA.to_numpy()
dma = eda.dudoit_mask(XAv)
La, ida = eda.to_log2(XAv)[:, dma], np.array(XA.columns)[dma]
syma = eda.probe_symbols("GPL8300", ida)
yb = binary_labels(y3, mll="ALL").to_numpy()


def collapse(L, syms, ids):
    """Her gen sembolu icin en yuksek varyansli probun sutun indeksi (etiket kullanmaz)."""
    var = L.var(axis=0)
    best = {}
    for j, s in enumerate(syms):
        if s == ids[j]:  # sembol bulunamadi
            continue
        if s not in best or var[j] > var[best[s]]:
            best[s] = j
    return best


bg, ba = collapse(Lg, symg, idg), collapse(La, syma, ida)
common = sorted(set(bg) & set(ba))
print(f"\nBirlesik analiz: Golub {len(bg)} sembol, Armstrong {len(ba)} sembol, ortak {len(common)} "
      f"(Dudoit filtresi sonrasi, prob basina ilk sembol)")
if len(common) < 50:
    sys.exit("Ortak gen cok az: once 09_common_genes.py calistirip platform tablolarini indir.")
G = Lg[:, [bg[s] for s in common]]
A = La[:, [ba[s] for s in common]]
dataset = np.array(["Golub"] * len(G) + ["Armstrong"] * len(A))
cls = np.concatenate([np.where(yg == 0, "ALL", "AML"), np.where(yb == 0, "ALL", "AML")])

variants = [
    ("Combined, global z-score", eda.zscore(np.vstack([G, A]))),
    ("Per-dataset z-score", np.vstack([eda.zscore(G), eda.zscore(A)])),
]
pal_d = {"Golub": "tab:green", "Armstrong": "tab:orange"}
pal_c = {"ALL": "tab:blue", "AML": "tab:red"}
fig, ax = plt.subplots(2, 2, figsize=(9, 7))
for r, (name, M) in enumerate(variants):
    scores, evr = eda.run_pca(M)
    sil_d = silhouette_score(scores[:, :2], dataset)
    sil_c = silhouette_score(scores[:, :2], cls)
    print(f"[{name}] PC1 = %{100 * evr[0]:.1f}, PC2 = %{100 * evr[1]:.1f}; "
          f"silhouette (veri seti) = {sil_d:.3f}, silhouette (sinif) = {sil_c:.3f}")
    for c, (lab_arr, pal, ttl) in enumerate([(dataset, pal_d, "dataset"), (cls, pal_c, "class")]):
        for lab, col in pal.items():
            m = lab_arr == lab
            ax[r, c].scatter(scores[m, 0], scores[m, 1], s=20, c=col, alpha=0.8, edgecolor="k",
                             linewidth=0.3, label=f"{lab} (n={m.sum()})")
        ax[r, c].set_title(f"{name} - {ttl}", fontsize=9)
        ax[r, c].set_xlabel(f"PC1 ({100 * evr[0]:.1f}%)")
        ax[r, c].set_ylabel(f"PC2 ({100 * evr[1]:.1f}%)")
        ax[r, c].legend(fontsize=7, frameon=False)
save(fig, "eda_joint_pca")
print("\nSekiller ve tablo ->", config.FIG_DIR)
print("Not: silhouette 1'e yakinsa gruplar ayrisiyor, 0'a yakinsa ic ice. Birlesik/global z-skorda veri seti "
      "silhouette'u yuksek, veri seti ici z-skorda dusukse veri seti farki ortalama duzeyden geliyor demektir.")

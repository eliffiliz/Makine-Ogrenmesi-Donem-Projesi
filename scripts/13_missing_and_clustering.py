"""EDA adim 5: (A) 'eksik/guvenilmez' deger ozeti (call sutunlari) ve (B) hiyerarsik kumeleme.

A) Golub ham CSV'lerindeki 'call' sutunlari (P=present, A=absent, M=marginal) incelenir:
   NaN/INF var mi, 'A' cagrilari ne kadar, dusuk degerler (<100) ile iliskisi ne.
B) Etiketsiz hiyerarsik kumeleme (Dudoit filtresi + log2 + gen z-skoru): farkli baglanti ve mesafeler,
   2 kume (Armstrong icin ayrica 3 kume) ile gercek siniflar arasindaki uyum (ARI, saflik).
   Etiketler yalnizca DEGERLENDIRMEDE kullanilir, kumeleme etiket gormez.
Ciktilar: figures/eda_golub_calls, eda_hclust_dendrogram (png+pdf), eda_golub_call_summary.csv,
          eda_hclust_results.csv
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_rand_score

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


# ================= A) call sutunlari ve eksik deger ozeti =================
def read_calls(path):
    raw = pd.read_csv(path)
    vals, calls, last = {}, {}, None
    for c in raw.columns:
        if c in ("Gene Description", "Gene Accession Number"):
            continue
        if str(c).startswith("call"):
            calls[last] = raw[c].astype(str).str.strip().to_numpy()
        else:
            last = int(c)
            vals[last] = pd.to_numeric(raw[c], errors="coerce").to_numpy()
    ids = raw["Gene Accession Number"].astype(str).to_numpy()
    return ids, pd.DataFrame(vals).to_numpy(), pd.DataFrame(calls).to_numpy()


rows, pooled_v, pooled_c = [], [], []
for name, path in [("train", config.TRAIN_FILE), ("test", config.TEST_FILE)]:
    ids, V, C = read_calls(path)
    keep = np.array([not i.startswith("AFFX-") for i in ids])
    V, C = V[keep], C[keep]
    low = V < 100
    absent = C == "A"
    row = {
        "split": name, "genes": V.shape[0], "samples": V.shape[1],
        "NaN": int(np.isnan(V).sum()), "INF": int(np.isinf(V).sum()),
        "P_ratio": float((C == "P").mean()), "A_ratio": float(absent.mean()), "M_ratio": float((C == "M").mean()),
        "A_and_value_below_100 / all_A": float((absent & low).sum() / max(absent.sum(), 1)),
        "value_below_100_and_A / all_below_100": float((absent & low).sum() / max(low.sum(), 1)),
        "genes_A_in_all_samples": float(absent.all(axis=1).mean()),
        "genes_A_in_90pct_samples": float((absent.mean(axis=1) >= 0.9).mean()),
    }
    rows.append(row)
    pooled_v.append(V.ravel())
    pooled_c.append(C.ravel())
summ = pd.DataFrame(rows)
summ.to_csv(config.FIG_DIR / "eda_golub_call_summary.csv", index=False)
print(summ.T.round(4).to_string(header=False))

pv, pc = np.concatenate(pooled_v), np.concatenate(pooled_c)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
bottoms = np.zeros(len(summ))
for cat, col in [("P_ratio", "tab:green"), ("M_ratio", "tab:orange"), ("A_ratio", "tab:red")]:
    ax[0].bar(summ["split"], summ[cat], bottom=bottoms, color=col, label={"P": "Present", "M": "Marginal", "A": "Absent"}[cat[0]])
    bottoms += summ[cat].to_numpy()
ax[0].set_ylabel("Fraction of measurements")
ax[0].set_title("Detection calls (P/M/A) by split", fontsize=9)
ax[0].legend(fontsize=7, loc="upper center", ncol=3, frameon=True, facecolor="white", framealpha=1)
L = np.log2(np.clip(pv, eda.FLOOR, eda.CEIL))
ax[1].boxplot([L[pc == c] for c in ["P", "M", "A"]], showfliers=False)
ax[1].set_xticklabels(["Present", "Marginal", "Absent"])
ax[1].set_ylabel(f"log2(clip[{eda.FLOOR:g}, {eda.CEIL:g}])")
ax[1].set_title("Expression value by detection call", fontsize=9)
save(fig, "eda_golub_calls")

# ================= B) hiyerarsik kumeleme =================
def prep(X_raw):
    dm = eda.dudoit_mask(X_raw)
    return eda.zscore(eda.to_log2(X_raw)[:, dm])


X_tr, y_tr, X_te, y_te = load_golub()
raw_ids = pd.read_csv(config.TRAIN_FILE, usecols=["Gene Accession Number"])["Gene Accession Number"].astype(str).to_numpy()
ctrl = np.array([i.startswith("AFFX-") for i in raw_ids])
Xg = pd.concat([X_tr, X_te]).to_numpy()[:, ~ctrl]
yg = pd.concat([y_tr, y_te]).to_numpy()
Mg = prep(Xg)

XA, y3 = load_armstrong()
Ma = prep(XA.to_numpy())
ya_bin = binary_labels(y3, mll="ALL").to_numpy()
ya3 = pd.factorize(y3.to_numpy())[0]

CONFIGS = [("ward", "euclidean"), ("average", "euclidean"), ("complete", "euclidean"),
           ("average", "correlation"), ("complete", "correlation")]


def purity(clusters, y):
    return sum(np.bincount(y[clusters == c]).max() for c in np.unique(clusters)) / len(y)


results = []
linkages = {}
for ds, M, tasks in [("Golub", Mg, [(2, yg, "ALL/AML")]),
                     ("Armstrong", Ma, [(2, ya_bin, "ALL(+MLL)/AML"), (3, ya3, "ALL/MLL/AML")])]:
    for method, metric in CONFIGS:
        Z = linkage(pdist(M, metric=metric), method=method)
        if method == "ward" and metric == "euclidean":
            linkages[ds] = Z
        for k, y, tname in tasks:
            cl = fcluster(Z, t=k, criterion="maxclust")
            results.append({"dataset": ds, "linkage": method, "metric": metric, "k": k, "task": tname,
                            "ARI": adjusted_rand_score(y, cl), "purity": purity(cl, y),
                            "cluster_sizes": sorted(np.bincount(cl)[1:].tolist(), reverse=True)})
res = pd.DataFrame(results)
res.to_csv(config.FIG_DIR / "eda_hclust_results.csv", index=False)
print("\nHiyerarsik kumeleme (etiketsiz; ARI 1 = tam uyum, 0 = rastgele):")
print(res.round(3).to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
for ax, (ds, y, pal, names) in zip(axes, [
        ("Golub", yg, {0: "tab:blue", 1: "tab:red"}, {0: "ALL", 1: "AML"}),
        ("Armstrong", ya3, {0: "tab:blue", 1: "tab:purple", 2: "tab:red"}, dict(enumerate(pd.factorize(y3.to_numpy())[1])))]):
    dn = dendrogram(linkages[ds], ax=ax, no_labels=True, color_threshold=0, above_threshold_color="gray")
    leaves = dn["leaves"]
    top = ax.get_ylim()[1]
    xs = 5 + 10 * np.arange(len(leaves))
    ax.scatter(xs, [-0.04 * top] * len(leaves), c=[pal[y[i]] for i in leaves], marker="s", s=14)
    ax.set_ylim(-0.08 * top, top)
    ax.set_title(f"{ds}: Ward linkage (unsupervised; colour = true class)", fontsize=9)
    ax.set_ylabel("Distance")
    for k, c in pal.items():
        ax.scatter([], [], c=c, marker="s", label=names[k])
    ax.legend(fontsize=7, frameon=False)
save(fig, "eda_hclust_dendrogram")
print("\nSekiller ve tablolar ->", config.FIG_DIR)

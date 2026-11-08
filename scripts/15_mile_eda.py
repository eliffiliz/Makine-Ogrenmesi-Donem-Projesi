"""EDA adim 7: MILE (GSE13159) ALL/AML alt kumesi icin EDA (Golub icin yaptigimizla ayni yapi).

Girdi (05_mile_expression.py uretir): data/mile/mile_expr_all_aml.npy, mile_probes.txt, mile_samples.csv
MILE degerleri ham Affymetrix sinyali DEGIL: GEO aciklamasina gore ozetleme + kuantil normalizasyon uygulanmis
ve degerler [0, 1] araliginda. Bu yuzden Dudoit esikleme/filtresi uygulanmaz; EDA donusumu log2(x + 0.01),
ardindan en yuksek varyansli 2000 prob ve z-skor (yalnizca kesif; etiket kullanmaz).
Ciktilar: figures/eda_mile_quality, eda_mile_projection, eda_mile_heatmap (png+pdf),
          mile_eda_summary.csv, mile_diff_probes.csv, mile_hclust_results.csv
DIKKAT: etiketli testler (Welch, isi haritasi) yalnizca kesiftir; modelleme/degerlendirmede kullanilmaz.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import mannwhitneyu
from sklearn.manifold import TSNE
from sklearn.metrics import adjusted_rand_score, silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config, eda

plt.rcParams.update({"font.size": 9})
config.FIG_DIR.mkdir(parents=True, exist_ok=True)
MILE = config.DATA_DIR / "mile"
OFFSET = 0.01
PAL_CLASS = {"ALL": "tab:blue", "AML": "tab:red"}
PAL_TISSUE = {"bone marrow": "tab:gray", "peripheral blood": "tab:olive"}


def save(fig, name):
    fig.savefig(config.FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(config.FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def scatter(ax, E, labels, palette, title, xl, yl, legend=True, small=False):
    labels = np.asarray(labels)
    for lab, col in palette.items():
        m = labels == lab
        if m.any():
            ax.scatter(E[m, 0], E[m, 1], s=8, c=[col], label=f"{lab} (n={m.sum()})", alpha=0.7, linewidth=0)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    if legend:
        if small:
            ax.legend(fontsize=5, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        else:
            ax.legend(fontsize=7, frameon=False)


X = np.load(MILE / "mile_expr_all_aml.npy")
probes = np.array((MILE / "mile_probes.txt").read_text(encoding="utf-8").split("\n"))
meta = pd.read_csv(MILE / "mile_samples.csv")
assert X.shape == (len(meta), len(probes)), "Matris, prob listesi ve ornek tablosu uyusmuyor"
y = meta["label"].to_numpy().astype(int)            # 0 = ALL, 1 = AML
cls = np.where(y == 0, "ALL", "AML")
tissue = meta["sample_type"].to_numpy()
subtype = meta["subtype"].to_numpy()
tissue_bin = (tissue == "peripheral blood").astype(int)
summary = {}

# ---------------- 1) veri kalitesi ----------------
nan_all = np.isnan(X).all(axis=0)
nan_part = np.isnan(X).any(axis=0) & ~nan_all
print(f"Ornek: {X.shape[0]}, prob: {X.shape[1]}; tamamen NaN prob: {int(nan_all.sum())}, kismen NaN prob: {int(nan_part.sum())}")
X, probes = X[:, ~nan_all], probes[~nan_all]
assert not np.isnan(X).any(), "Kismen NaN iceren problar var; ayri islem gerekir"
rng = np.random.default_rng(config.SEED)
v = rng.choice(X.ravel(), size=min(2_000_000, X.size), replace=False)
pct = np.percentile(v, [1, 50, 99])
summary.update({"samples": X.shape[0], "probes_after_NaN_removal": X.shape[1], "all_NaN_probes": int(nan_all.sum()),
                "min": float(X.min()), "max": float(X.max()), "p1": float(pct[0]), "median": float(pct[1]),
                "p99": float(pct[2]), "fraction_zero": float((v == 0).mean()),
                "median_of_sample_medians": float(np.median(np.median(X, axis=1))),
                "sd_of_sample_medians": float(np.median(X, axis=1).std())})
print(f"Deger araligi [{X.min():.3f}, {X.max():.3f}], medyan {pct[1]:.3f}; sifir orani {(v == 0).mean():.4f}; "
      f"ornek medyanlarinin std'si {np.median(X, axis=1).std():.2e} (kuantil normalizasyon: ~0 beklenir)")
ct = pd.crosstab(cls, tissue)
print("Sinif x doku (adet):\n", ct.to_string())

fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))
bottom = np.zeros(2)
for t, col in PAL_TISSUE.items():
    vals = np.array([int(((cls == c) & (tissue == t)).sum()) for c in ["ALL", "AML"]])
    ax[0].bar(["ALL", "AML"], vals, bottom=bottom, color=col, label=t)
    bottom += vals
ax[0].set_ylabel("Number of samples")
ax[0].set_title("Class and sample source", fontsize=9)
ax[0].legend(fontsize=7, frameon=False)
ax[1].hist(v, bins=100, density=True)
ax[1].set_title("Values as provided (GEO matrix)", fontsize=9)
ax[1].set_xlabel("Expression value")
ax[1].set_ylabel("Density")
ax[2].hist(np.log2(v + OFFSET), bins=100, density=True)
ax[2].set_title(f"log2(x + {OFFSET:g})", fontsize=9)
ax[2].set_xlabel("log2 expression")
fig.tight_layout()
save(fig, "eda_mile_quality")

# ---------------- 2) donusum ve ozellik secimi ----------------
L = np.log2(X + OFFSET).astype(np.float32)
idx = eda.top_variance(L, 2000)
Z = eda.zscore(L[:, idx].astype(np.float64))
scores, evr = eda.run_pca(Z)
print(f"PCA (en yuksek varyansli {len(idx)} prob): PC1 = %{100 * evr[0]:.1f}, PC2 = %{100 * evr[1]:.1f}, ilk 10 PC = %{100 * evr.sum():.1f}")
tsne = TSNE(n_components=2, perplexity=30, init="pca", random_state=config.SEED).fit_transform(Z)
try:
    import umap
    um = umap.UMAP(n_neighbors=15, min_dist=0.3, random_state=config.SEED).fit_transform(Z)
except Exception as e:
    print("UMAP atlandi:", e)
    um = None

pal_sub = {s: plt.cm.tab20(i % 20) for i, s in enumerate(sorted(set(subtype)))}
fig, ax = plt.subplots(2, 3, figsize=(12, 7))
pcx, pcy = f"PC1 ({100 * evr[0]:.1f}%)", f"PC2 ({100 * evr[1]:.1f}%)"
scatter(ax[0, 0], scores[:, :2], cls, PAL_CLASS, "PCA - class", pcx, pcy)
scatter(ax[0, 1], scores[:, :2], tissue, PAL_TISSUE, "PCA - sample source", pcx, pcy)
ax[0, 2].bar(range(1, len(evr) + 1), 100 * evr)
ax[0, 2].set_title("Explained variance", fontsize=9)
ax[0, 2].set_xlabel("Principal component")
ax[0, 2].set_ylabel("Variance (%)")
scatter(ax[1, 0], tsne, cls, PAL_CLASS, "t-SNE - class", "t-SNE 1", "t-SNE 2")
if um is not None:
    scatter(ax[1, 1], um, cls, PAL_CLASS, "UMAP - class", "UMAP 1", "UMAP 2")
    scatter(ax[1, 2], um, subtype, pal_sub, "UMAP - subtype", "UMAP 1", "UMAP 2", small=True)
else:
    ax[1, 1].axis("off")
    ax[1, 2].axis("off")
fig.tight_layout()
save(fig, "eda_mile_projection")

# ---------------- 3) ayirt edici problar (kesif) ----------------
syms = np.array(eda.probe_symbols("GPL570", probes))
de = eda.diff_expression(L.astype(np.float64), y, probes, syms)
de["abs_t"] = de["t"].abs()
de = de.sort_values("abs_t", ascending=False).reset_index(drop=True)
sig = de["fdr"] < 0.05
up = float((de.loc[sig, "log2fc"] > 0).mean())
summary.update({"probes_tested": len(de), "fdr_lt_0.05": int(sig.sum()), "fdr_lt_0.05_fraction": float(sig.mean()),
                "fdr_lt_0.05_and_abs_log2fc_gt_1": int((sig & (de["log2fc"].abs() > 1)).sum()),
                "sig_higher_in_AML_fraction": up})
print(f"Welch: {len(de)} prob test edildi, FDR<0.05: {int(sig.sum())} (%{100 * sig.mean():.1f}); "
      f"|log2FC|>1: {int((sig & (de['log2fc'].abs() > 1)).sum())}; anlamli problarin %{100 * up:.1f}'i AML'de yuksek")
print("|t| en buyuk 5 prob (n=1292 oldugu icin p-degerleri 0'a yakin, siralama |t| ile):")
print(de.head(5)[["probe", "symbol", "t", "log2fc"]].round(3).to_string(index=False))
de.drop(columns=["col_index"]).head(100).to_csv(config.FIG_DIR / "mile_diff_probes.csv", index=False)

top = de.head(40)
sel = np.concatenate([rng.choice(np.where(y == 0)[0], 60, replace=False),
                      rng.choice(np.where(y == 1)[0], 60, replace=False)])
Zh = eda.zscore(L[sel][:, top["col_index"].to_numpy()].astype(np.float64))
rows = [f"{s} ({p})" for s, p in zip(top["symbol"], top["probe"])]
sid = [f"S{i}" for i in sel]
data = pd.DataFrame(Zh.T, index=rows, columns=sid)
colors = pd.DataFrame({"class": [PAL_CLASS[c] for c in cls[sel]],
                       "source": [PAL_TISSUE[t] for t in tissue[sel]]}, index=sid)
cg = sns.clustermap(data, method="average", metric="euclidean", cmap="vlag", center=0, vmin=-3, vmax=3,
                    col_colors=colors, figsize=(10, 9), xticklabels=False, yticklabels=True,
                    cbar_kws={"label": "z-score"})
cg.ax_heatmap.tick_params(axis="y", labelsize=6)
cg.fig.suptitle("MILE: top 40 probes by |t| (exploratory, labels used; 60 ALL + 60 AML random samples)", y=1.01, fontsize=9)
cg.fig.savefig(config.FIG_DIR / "eda_mile_heatmap.png", dpi=300, bbox_inches="tight")
cg.fig.savefig(config.FIG_DIR / "eda_mile_heatmap.pdf", bbox_inches="tight")
plt.close(cg.fig)

# ---------------- 4) etiketsiz hiyerarsik kumeleme ----------------
def purity(cl, labels):
    return sum(np.bincount(labels[cl == c]).max() for c in np.unique(cl)) / len(labels)

res = []
for method, metric in [("ward", "euclidean"), ("average", "euclidean"), ("complete", "euclidean"),
                       ("average", "correlation"), ("complete", "correlation")]:
    cl = fcluster(linkage(pdist(Z, metric=metric), method=method), t=2, criterion="maxclust")
    res.append({"linkage": method, "metric": metric, "ARI_vs_class": adjusted_rand_score(y, cl),
                "purity_vs_class": purity(cl, y), "ARI_vs_source": adjusted_rand_score(tissue_bin, cl),
                "cluster_sizes": sorted(np.bincount(cl)[1:].tolist(), reverse=True)})
res = pd.DataFrame(res)
res.to_csv(config.FIG_DIR / "mile_hclust_results.csv", index=False)
print("\nEtiketsiz hiyerarsik kumeleme, k=2 (ARI: 1 tam uyum, 0 rastgele):")
print(res.round(3).to_string(index=False))

# ---------------- 5) doku kaynagi (batch ipucu) ----------------
s_cls, s_tis = silhouette_score(scores[:, :2], y), silhouette_score(scores[:, :2], tissue_bin)
summary.update({"silhouette_PC12_by_class": float(s_cls), "silhouette_PC12_by_source": float(s_tis)})
print(f"\nSilhouette (PC1-PC2): sinifa gore {s_cls:.3f}, doku kaynagina gore {s_tis:.3f}")
for name, c in [("ALL", 0), ("AML", 1)]:
    m = y == c
    for k in (0, 1):
        p = mannwhitneyu(scores[m & (tissue_bin == 0), k], scores[m & (tissue_bin == 1), k]).pvalue
        print(f"[{name}] PC{k + 1}: kemik iligi (n={int((m & (tissue_bin == 0)).sum())}) vs kan (n={int((m & (tissue_bin == 1)).sum())}), Mann-Whitney p = {p:.3g}")
    dt = eda.diff_expression(L[m].astype(np.float64), tissue_bin[m], probes, syms)
    frac = float((dt["fdr"] < 0.05).mean())
    summary[f"{name}_BM_vs_PB_sig_fraction"] = frac
    print(f"[{name}] kemik iligi vs kan: {len(dt)} prob test edildi, FDR<0.05 orani = %{100 * frac:.1f} "
          f"(sinif farki icin: %{100 * sig.mean():.1f})")

pd.DataFrame([summary]).T.rename(columns={0: "value"}).to_csv(config.FIG_DIR / "mile_eda_summary.csv")
print("\nSekiller ve tablolar ->", config.FIG_DIR)

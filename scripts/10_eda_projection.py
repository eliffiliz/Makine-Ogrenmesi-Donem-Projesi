"""EDA adim 2: PCA, t-SNE, UMAP, ayirt edici genler, isi haritasi + hiyerarsik kumeleme.

Golub (72 ornek) ve Armstrong (72 ornek) icin ayri ayri calisir.
Onisleme (yalnizca EDA): sabit esiklerle kirp [100, 16000] -> log2 -> en yuksek varyansli
2000 gen -> z-skor. Etiketler yalniz renklendirme ve fark testi icin kullanilir.

DIKKAT: Etiketle secilen genler (fark testi, isi haritasi) yalnizca gorsellestirme icindir;
modelleme/degerlendirmede kullanilmaz.
Cikti: figures/eda_<veri>_projection, eda_<veri>_heatmap (png+pdf), <veri>_diff_genes.csv
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from sklearn.manifold import TSNE

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config, eda
from src.armstrong_loading import binary_labels, load_armstrong
from src.data_loading import load_golub

plt.rcParams.update({"font.size": 9})
config.FIG_DIR.mkdir(parents=True, exist_ok=True)
PAL_CLASS = {"ALL": "tab:blue", "AML": "tab:red"}
PAL_SPLIT = {"train": "tab:green", "test": "tab:orange"}
PAL_SUB = {"ALL": "tab:blue", "MLL": "tab:purple", "AML": "tab:red"}


def save(fig, name):
    fig.savefig(config.FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(config.FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def scatter(ax, E, labels, palette, title, xl, yl):
    labels = np.asarray(labels)
    for lab, col in palette.items():
        m = labels == lab
        ax.scatter(E[m, 0], E[m, 1], s=24, c=col, label=f"{lab} (n={m.sum()})",
                   alpha=0.85, edgecolor="k", linewidth=0.3)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.legend(fontsize=7, frameon=False)


def umap_embed(M):
    try:
        import umap
        return umap.UMAP(n_neighbors=10, min_dist=0.3, random_state=config.SEED).fit_transform(M)
    except Exception as e:  # umap kurulu degilse ya da hata verirse atla
        print("UMAP atlandi:", e)
        return None


def analyze(name, L, y_bin, class_str, aux, pal_aux, aux_title, probe_ids, symbols, sample_ids):
    print(f"\n===== {name}: {L.shape[0]} ornek x {L.shape[1]} gen =====")
    idx = eda.top_variance(L, 2000)
    M = eda.zscore(L[:, idx])
    scores, evr = eda.run_pca(M)
    print(f"PCA: PC1 = %{100 * evr[0]:.1f}, PC2 = %{100 * evr[1]:.1f}, ilk 10 PC toplami = %{100 * evr.sum():.1f}")
    tsne = TSNE(n_components=2, perplexity=15, init="pca", random_state=config.SEED).fit_transform(M)
    um = umap_embed(M)

    fig, ax = plt.subplots(2, 3, figsize=(11, 6.5))
    pcx, pcy = f"PC1 ({100 * evr[0]:.1f}%)", f"PC2 ({100 * evr[1]:.1f}%)"
    scatter(ax[0, 0], scores[:, :2], class_str, PAL_CLASS, "PCA - class", pcx, pcy)
    scatter(ax[0, 1], scores[:, :2], aux, pal_aux, f"PCA - {aux_title}", pcx, pcy)
    ax[0, 2].bar(range(1, len(evr) + 1), 100 * evr)
    ax[0, 2].set_title("Explained variance", fontsize=9)
    ax[0, 2].set_xlabel("Principal component")
    ax[0, 2].set_ylabel("Variance (%)")
    scatter(ax[1, 0], tsne, class_str, PAL_CLASS, "t-SNE - class", "t-SNE 1", "t-SNE 2")
    if um is not None:
        scatter(ax[1, 1], um, class_str, PAL_CLASS, "UMAP - class", "UMAP 1", "UMAP 2")
    else:
        ax[1, 1].axis("off")
    scatter(ax[1, 2], tsne, aux, pal_aux, f"t-SNE - {aux_title}", "t-SNE 1", "t-SNE 2")
    fig.tight_layout()
    save(fig, f"eda_{name}_projection")

    # --- ayirt edici genler (yalnizca kesif) ---
    de = eda.diff_expression(L, y_bin, probe_ids, symbols)
    sig = de["fdr"] < 0.05
    sig_fc = sig & (de["log2fc"].abs() > 1)
    print(f"Test edilen gen: {len(de)}, FDR<0.05: {int(sig.sum())}, FDR<0.05 ve |log2FC|>1: {int(sig_fc.sum())}")
    print("En iyi 5 gen:")
    print(de.head(5)[["probe", "symbol", "t", "fdr", "log2fc"]].round(4).to_string(index=False))
    de.drop(columns="col_index").head(100).to_csv(config.FIG_DIR / f"{name}_diff_genes.csv", index=False)

    top = de.head(40)
    Z = eda.zscore(L[:, top["col_index"].to_numpy()])
    rows = [f"{s} ({p})" for s, p in zip(top["symbol"], top["probe"])]
    data = pd.DataFrame(Z.T, index=rows, columns=sample_ids)
    cols = pd.Series([PAL_CLASS[c] for c in class_str], index=sample_ids)
    cg = sns.clustermap(data, method="average", metric="euclidean", cmap="vlag", center=0,
                        vmin=-3, vmax=3, col_colors=cols, figsize=(10, 9),
                        xticklabels=False, yticklabels=True, cbar_kws={"label": "z-score"})
    cg.ax_heatmap.tick_params(axis="y", labelsize=6)
    cg.fig.suptitle(f"{name}: top 40 genes by Welch t-test (exploratory, labels used)", y=1.01, fontsize=9)
    cg.fig.savefig(config.FIG_DIR / f"eda_{name}_heatmap.png", dpi=300, bbox_inches="tight")
    cg.fig.savefig(config.FIG_DIR / f"eda_{name}_heatmap.pdf", bbox_inches="tight")
    plt.close(cg.fig)
    return scores, evr


# ---------------- Golub ----------------
X_tr, y_tr, X_te, y_te = load_golub()
raw_ids = pd.read_csv(config.TRAIN_FILE, usecols=["Gene Accession Number"])["Gene Accession Number"].astype(str).to_numpy()
assert len(raw_ids) == X_tr.shape[1]
keep = np.array([not i.startswith("AFFX-") for i in raw_ids])
X_all = pd.concat([X_tr, X_te])
y_all = pd.concat([y_tr, y_te])
Xk = X_all.to_numpy()[:, keep]
dm = eda.dudoit_mask(Xk)
print(f"Golub: {Xk.shape[1]} prob -> Dudoit filtresi sonrasi {int(dm.sum())}")
L = eda.to_log2(Xk)[:, dm]
ids = raw_ids[keep][dm]
syms = np.array(eda.probe_symbols("GPL80", ids))
split = np.array(["train"] * len(X_tr) + ["test"] * len(X_te))
cls = np.where(y_all.to_numpy() == 0, "ALL", "AML")
scores, evr = analyze("golub", L, y_all.to_numpy(), cls, split, PAL_SPLIT, "train/test split",
                      ids, syms, [f"P{i}" for i in X_all.index])
# batch ipucu: egitim ve test arasinda PC skorlari ve ornek medyani
med = np.median(L, axis=1)
for label, v in [("PC1", scores[:, 0]), ("PC2", scores[:, 1]), ("ornek medyani (log2)", med)]:
    p = mannwhitneyu(v[split == "train"], v[split == "test"]).pvalue
    print(f"Egitim vs test, {label}: Mann-Whitney p = {p:.3g}")
print("Not: egitim/test sinif orani farkli (egitim %71 ALL, test %59 ALL); sinif etkisi karisabilir.")

# ---------------- Armstrong ----------------
XA, y3 = load_armstrong()
yb = binary_labels(y3, mll="ALL").to_numpy()
XAv = XA.to_numpy()
dmA = eda.dudoit_mask(XAv)
print(f"Armstrong: {XAv.shape[1]} prob -> Dudoit filtresi sonrasi {int(dmA.sum())}")
LA = eda.to_log2(XAv)[:, dmA]
idsA = np.array(XA.columns)[dmA]
symsA = np.array(eda.probe_symbols("GPL8300", idsA))
clsA = np.where(yb == 0, "ALL", "AML")
analyze("armstrong", LA, yb, clsA, y3.to_numpy(), PAL_SUB, "subtype",
        idsA, symsA, list(XA.index))
print("\nSekiller ve tablolar ->", config.FIG_DIR)

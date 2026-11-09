"""EDA adim 6: Golub egitim ve bagimsiz test kumeleri birbirinden ayirt edilebiliyor mu?

Raporun 'egitim/test ayrimi' paragrafindaki TUM sayilari uretir:
 1) ornek medyani farki (test - egitim): tum ornekler, yalniz ALL, yalniz AML (Mann-Whitney)
 2) egitim/test etiketini tahmin eden L2-lojistik regresyonun tekrarli tabakali 5 katli ROC-AUC'si
    (olcekleme kat icinde ogrenilir; ozellik: Dudoit filtresi sonrasi log2 genler)
 3) permutasyon testi: egitim/test etiketleri SINIF ICINDE karistirilir (sinif bilesimi sabit kalir)
 4) sinif icinde egitim vs test gen testi (Welch + BH-FDR): anlamli gen orani
 5) PCA'nin ilk iki bileseninde silhouette (egitim/test'e gore ve sinifa gore)
Ciktilar: figures/eda_split_distinguishability.csv, figures/eda_split_null_auc.(png|pdf)
"""
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, silhouette_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config, eda
from src.data_loading import load_golub

warnings.filterwarnings("ignore")
N_PERM = 200
SEED = config.SEED
config.FIG_DIR.mkdir(parents=True, exist_ok=True)

X_tr, y_tr, X_te, y_te = load_golub()
raw_ids = pd.read_csv(config.TRAIN_FILE, usecols=["Gene Accession Number"])["Gene Accession Number"].astype(str).to_numpy()
ctrl = np.array([i.startswith("AFFX-") for i in raw_ids])
Xg = pd.concat([X_tr, X_te]).to_numpy()[:, ~ctrl]
cls = pd.concat([y_tr, y_te]).to_numpy()                      # 0 = ALL, 1 = AML
split = np.array([0] * len(X_tr) + [1] * len(X_te))           # 0 = egitim, 1 = test
dm = eda.dudoit_mask(Xg)
L = eda.to_log2(Xg)[:, dm]
ids = raw_ids[~ctrl][dm]
print(f"Ozellik: {L.shape[1]} gen (Dudoit filtresi sonrasi), {L.shape[0]} ornek; N_PERM = {N_PERM}")

rows = []

# ---- 1) ornek medyani farki ----
med = np.median(L, axis=1)
for name, m in [("all", np.ones(len(med), bool)), ("ALL only", cls == 0), ("AML only", cls == 1)]:
    a, b = med[m & (split == 0)], med[m & (split == 1)]
    p = mannwhitneyu(a, b).pvalue
    rows.append({"analysis": "median shift (test - train)", "subset": name,
                 "n_train": len(a), "n_test": len(b),
                 "mean_diff": float(b.mean() - a.mean()), "median_diff": float(np.median(b) - np.median(a)),
                 "p_value": float(p)})
    print(f"[ornek medyani] {name}: test-egitim ortalama fark = {b.mean() - a.mean():+.3f}, "
          f"medyan fark = {np.median(b) - np.median(a):+.3f}, Mann-Whitney p = {p:.3g}")


# ---- 2-3) siniflandirici AUC + sinif ici permutasyon ----
def cv_auc(M, y, repeats, seed):
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=repeats, random_state=seed)
    out = []
    for tr, te in rskf.split(M, y):
        clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000))
        clf.fit(M[tr], y[tr])
        out.append(roc_auc_score(y[te], clf.predict_proba(M[te])[:, 1]))
    return np.array(out)


def permute_within_class(y_split, y_class, rng):
    out = y_split.copy()
    for c in np.unique(y_class):
        idx = np.where(y_class == c)[0]
        out[idx] = rng.permutation(out[idx])
    return out


nulls = {}
subsets = [("all genes, all samples", np.ones(len(cls), bool), None),
           ("all genes, ALL only", cls == 0, None),
           ("all genes, AML only", cls == 1, None),
           ("top-2000 variance genes, all samples", np.ones(len(cls), bool), 2000)]
for name, m, topk in subsets:
    M = L[m]
    if topk:
        M = M[:, eda.top_variance(M, topk)]
    y, c = split[m], cls[m]
    obs = cv_auc(M, y, repeats=5, seed=SEED)                       # raporlanan: 5 tekrar x 5 kat
    stat = cv_auc(M, y, repeats=1, seed=SEED).mean()               # permutasyon istatistigi (1 tekrar)
    rng = np.random.default_rng(SEED)
    null = np.array([cv_auc(M, permute_within_class(y, c, rng), repeats=1, seed=SEED).mean() for _ in range(N_PERM)])
    p = (1 + int((null >= stat).sum())) / (1 + N_PERM)
    nulls[name] = (null, stat)
    rows.append({"analysis": "train-vs-test classifier AUC", "subset": name, "n_train": int((y == 0).sum()),
                 "n_test": int((y == 1).sum()), "auc_mean": float(obs.mean()), "auc_std": float(obs.std()),
                 "null_mean": float(null.mean()), "perm_p": float(p)})
    print(f"[AUC] {name}: {obs.mean():.3f} +- {obs.std():.3f} (5x5 kat); permutasyon bos dagilim ort. = "
          f"{null.mean():.3f}; permutasyon p = {p:.3f}")

# ---- 4) sinif ici gen testi ----
for name, c in [("ALL only", 0), ("AML only", 1)]:
    m = cls == c
    de = eda.diff_expression(L[m], split[m], ids, ids)
    frac = float((de["fdr"] < 0.05).mean())
    rows.append({"analysis": "within-class train-vs-test gene test", "subset": name,
                 "n_train": int((split[m] == 0).sum()), "n_test": int((split[m] == 1).sum()),
                 "genes_tested": len(de), "sig_fraction": frac})
    print(f"[gen testi] {name}: {len(de)} gen test edildi, FDR<0.05 orani = %{100 * frac:.2f}")

# ---- 5) PCA'da silhouette ----
M = eda.zscore(L)
scores, evr = eda.run_pca(M)
s_split = silhouette_score(scores[:, :2], split)
s_cls = silhouette_score(scores[:, :2], cls)
rows.append({"analysis": "PCA silhouette", "subset": "all", "silhouette_by_split": float(s_split),
             "silhouette_by_class": float(s_cls)})
print(f"[silhouette, PC1-PC2] egitim/test'e gore = {s_split:.3f}, sinifa gore = {s_cls:.3f}")

pd.DataFrame(rows).to_csv(config.FIG_DIR / "eda_split_distinguishability.csv", index=False)

fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
for ax, name in zip(axes, [s[0] for s in subsets[:3]]):
    null, stat = nulls[name]
    ax.hist(null, bins=25, color="lightgray", edgecolor="k")
    ax.axvline(stat, color="tab:red", label=f"observed (1x5-fold CV) = {stat:.2f}")
    ax.set_title(name, fontsize=9)
    ax.set_xlabel("Mean CV ROC-AUC (permuted labels)")
    ax.set_ylabel("Count")
    ax.legend(fontsize=7, frameon=False)
fig.tight_layout()
fig.savefig(config.FIG_DIR / "eda_split_null_auc.png", dpi=300)
fig.savefig(config.FIG_DIR / "eda_split_null_auc.pdf")
print("\nTablo ve sekil ->", config.FIG_DIR)

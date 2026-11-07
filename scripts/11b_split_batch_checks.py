"""EDA adim 3b: Golub egitim (n=38) ve bagimsiz test (n=34) arasi ayirt edilebilirlik (batch) kontrolu.

Soru: egitim ve test kumeleri ifade duzeyinde birbirinden ayirt edilebiliyor mu?
Sinif bilesimi iki kumede farkli oldugu icin (egitim 27 ALL / 11 AML, test 20 ALL / 14 AML)
fark, sinif icinde de (yalnizca ALL, yalnizca AML) ayri ayri olculur.

Yontemler (hepsi kesif amaclidir, model egitimi/degerlendirmesinde KULLANILMAZ):
  1) Ornek basina medyan log2 ifade: egitim vs test (Mann-Whitney U), sinif icinde.
  2) Gen basina Welch t-testi egitim vs test (BH-FDR < 0.05 olan gen orani), sinif icinde.
  3) "Kume siniflandirici": egitim/test etiketini tahmin eden L2-lojistik regresyonun
     tekrarli tabakali CV ROC-AUC'si. ~0.5 ayirt edilemiyor, 1'e yakin belirgin batch.
  4) PCA: silhouette (kume ile) ve silhouette (sinif ile).

Not: gen filtresi (Dudoit) ve en yuksek varyansli 2000 gen secimi ne sinif ne de kume
etiketi kullanir; yine de bu sayilar yalnizca EDA icindir.
Dikkat: yalnizca AML icinde kumeler cok kucuk (11 vs 14), bu yuzden CV AUC'si cok oynak
olur ve 0.5'in altina da dusebilir; yorumlarken std degerine bakin, tek basina AUC'ye degil.

Cikti: figures/eda_golub_batch.(png|pdf), figures/eda_golub_batch_summary.csv
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import silhouette_score
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config, eda
from src.data_loading import load_golub

plt.rcParams.update({"font.size": 9})
config.FIG_DIR.mkdir(parents=True, exist_ok=True)
N_TOP_GENES = 10**6
FDR_ALPHA = 0.05


def save(fig, name):
    fig.tight_layout()
    fig.savefig(config.FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(config.FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------- veri ----------------
X_tr, y_tr, X_te, y_te = load_golub()
raw_ids = pd.read_csv(config.TRAIN_FILE, usecols=["Gene Accession Number"])[
    "Gene Accession Number"].astype(str).to_numpy()
ctrl = np.array([i.startswith("AFFX-") for i in raw_ids])  # Affymetrix kontrol problari

X = np.vstack([X_tr.to_numpy(), X_te.to_numpy()])[:, ~ctrl]
cls = np.concatenate([y_tr.to_numpy(), y_te.to_numpy()])          # 0 = ALL, 1 = AML
split = np.r_[np.zeros(len(X_tr), int), np.ones(len(X_te), int)]  # 0 = egitim, 1 = test

keep = eda.dudoit_mask(X)
L = eda.to_log2(X)[:, keep]
print(f"Ornek: {L.shape[0]} (egitim {int((split == 0).sum())}, test {int((split == 1).sum())}); "
      f"gen: {L.shape[1]} (kontrol probsuz {X.shape[1]}, Dudoit filtresi sonrasi)")

top = eda.top_variance(L, k=N_TOP_GENES)
Ltop = L[:, top]
scopes = {"Tum ornekler": np.ones(len(cls), bool),
          "Yalniz ALL": cls == 0,
          "Yalniz AML": cls == 1}

# ---------------- 1-3) sinif icinde / tumde karsilastirma ----------------
med = np.median(L, axis=1)
rows = []
auc_mean, auc_std = {}, {}
for name, m in scopes.items():
    s = split[m]
    n_tr, n_te = int((s == 0).sum()), int((s == 1).sum())

    u = stats.mannwhitneyu(med[m][s == 0], med[m][s == 1], alternative="two-sided")
    med_diff = float(np.median(med[m][s == 1]) - np.median(med[m][s == 0]))

    t, p = stats.ttest_ind(L[m][s == 1], L[m][s == 0], equal_var=False)
    ok = ~np.isnan(p)
    fdr = eda.bh_fdr(p[ok])
    frac_sig = float((fdr < FDR_ALPHA).mean())

    k = min(5, n_tr, n_te)
    cv = RepeatedStratifiedKFold(n_splits=k, n_repeats=5, random_state=config.SEED)
    pipe = make_pipeline(StandardScaler(),
                         LogisticRegression(C=0.1, max_iter=2000, random_state=config.SEED))
    aucs = cross_val_score(pipe, Ltop[m], s, cv=cv, scoring="roc_auc")
    auc_mean[name], auc_std[name] = float(aucs.mean()), float(aucs.std())

    rows.append({"kapsam": name, "n_egitim": n_tr, "n_test": n_te,
                 "medyan_fark_test_eksi_egitim": med_diff, "mannwhitney_p": float(u.pvalue),
                 f"gen_orani_FDR<{FDR_ALPHA}": frac_sig,
                 "kume_siniflandirici_AUC_ort": auc_mean[name],
                 "kume_siniflandirici_AUC_std": auc_std[name], "cv_kat": k})
    print(f"[{name}] egitim={n_tr}, test={n_te} | medyan fark={med_diff:+.3f} (MWU p={u.pvalue:.3g}) | "
          f"FDR<{FDR_ALPHA} gen orani={frac_sig:.3f} | kume AUC={aucs.mean():.3f} +- {aucs.std():.3f}")

# ---------------- 4) PCA ----------------
scores, evr = eda.run_pca(eda.zscore(Ltop))
sil_split = silhouette_score(scores[:, :2], split)
sil_cls = silhouette_score(scores[:, :2], cls)
print(f"PCA: PC1 = %{100 * evr[0]:.1f}, PC2 = %{100 * evr[1]:.1f}; "
      f"silhouette (egitim/test) = {sil_split:.3f}, silhouette (sinif) = {sil_cls:.3f}")
rows.append({"kapsam": "PCA (ilk 2 PC, en yuksek varyansli 2000 gen)",
             "silhouette_kume": sil_split, "silhouette_sinif": sil_cls,
             "PC1_varyans": float(evr[0]), "PC2_varyans": float(evr[1])})
pd.DataFrame(rows).to_csv(config.FIG_DIR / "eda_golub_batch_summary.csv", index=False)

# ---------------- sekil ----------------
fig, ax = plt.subplots(2, 2, figsize=(9, 7))
cs = {0: "tab:blue", 1: "tab:orange"}
for v, lab in [(0, "Train (n=38)"), (1, "Independent test (n=34)")]:
    mm = split == v
    ax[0, 0].scatter(scores[mm, 0], scores[mm, 1], s=22, c=cs[v], alpha=0.8,
                     edgecolor="k", linewidth=0.3, label=lab)
ax[0, 0].set_title(f"PCA by split (silhouette = {sil_split:.2f})", fontsize=9)
cc = {0: "tab:blue", 1: "tab:red"}
for v, lab in [(0, "ALL"), (1, "AML")]:
    mm = cls == v
    ax[0, 1].scatter(scores[mm, 0], scores[mm, 1], s=22, c=cc[v], alpha=0.8,
                     edgecolor="k", linewidth=0.3, label=lab)
ax[0, 1].set_title(f"PCA by class (silhouette = {sil_cls:.2f})", fontsize=9)
for a in ax[0]:
    a.set_xlabel(f"PC1 ({100 * evr[0]:.1f}%)")
    a.set_ylabel(f"PC2 ({100 * evr[1]:.1f}%)")
    a.legend(fontsize=7, frameon=False)

groups, names = [], []
for c, cn in [(0, "ALL"), (1, "AML")]:
    for v, vn in [(0, "train"), (1, "test")]:
        groups.append(med[(cls == c) & (split == v)])
        names.append(f"{cn}\n{vn}")
ax[1, 0].boxplot(groups, showmeans=True)
ax[1, 0].set_xticklabels(names)
ax[1, 0].set_ylabel("Per-sample median log2 expression")
ax[1, 0].set_title("Sample medians by class and split", fontsize=9)

labs = list(scopes)
ax[1, 1].bar(range(len(labs)), [auc_mean[k] for k in labs],
             yerr=[auc_std[k] for k in labs], capsize=4, color="tab:gray")
ax[1, 1].axhline(0.5, color="k", linestyle="--", linewidth=0.8)
ax[1, 1].set_xticks(range(len(labs)))
ax[1, 1].set_xticklabels(["All", "ALL only", "AML only"])
ax[1, 1].set_ylim(0, 1.05)
ax[1, 1].set_ylabel("CV ROC-AUC (predict train vs test)")
ax[1, 1].set_title("Can a classifier tell the splits apart?", fontsize=9)
save(fig, "eda_golub_batch")
print("Sekil ve tablo ->", config.FIG_DIR)
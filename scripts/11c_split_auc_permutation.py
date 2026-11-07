"""EDA adim 3c: egitim/test "kume siniflandirici" AUC'si icin permutasyon bos dagilimi.

11b'de egitim/test etiketini tahmin eden siniflandiricinin CV ROC-AUC'si ~0.77-0.81 cikti.
Bu betik o degerin TESADUFEN ne kadar olasi oldugunu olcer:
  - egitim/test etiketleri SINIF ICINDE rastgele karistirilir (her sinifin egitim/test sayisi
    ve dolayisiyla sinif bilesimi aynen korunur),
  - her permutasyonda ayni siniflandirici ayni CV ile degerlendirilir,
  - tek tarafli p = (1 + #{bos AUC >= gozlenen AUC}) / (1 + n_perm).

Kapsamlar: "Tum ornekler" (sinif bilesimi sabit tutularak karistirilir), "Yalniz ALL", "Yalniz AML".
Yalnizca kesif amaclidir; model egitimi/degerlendirmesinde kullanilmaz.

Cikti: figures/eda_golub_split_auc_permutation.csv (+ konsol)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config, eda
from src.data_loading import load_golub

N_PERM = 200
N_TOP_GENES = 10**6

X_tr, y_tr, X_te, y_te = load_golub()
raw_ids = pd.read_csv(config.TRAIN_FILE, usecols=["Gene Accession Number"])[
    "Gene Accession Number"].astype(str).to_numpy()
ctrl = np.array([i.startswith("AFFX-") for i in raw_ids])

X = np.vstack([X_tr.to_numpy(), X_te.to_numpy()])[:, ~ctrl]
cls = np.concatenate([y_tr.to_numpy(), y_te.to_numpy()])
split = np.r_[np.zeros(len(X_tr), int), np.ones(len(X_te), int)]

L = eda.to_log2(X)[:, eda.dudoit_mask(X)]
Ltop = L[:, eda.top_variance(L, k=N_TOP_GENES)]
scopes = {"Tum ornekler": np.ones(len(cls), bool), "Yalniz ALL": cls == 0, "Yalniz AML": cls == 1}


def make_pipe():
    return make_pipeline(StandardScaler(),
                         LogisticRegression(C=0.1, max_iter=2000, random_state=config.SEED))


def cv_auc(Z, s, k, repeats, seed):
    cv = RepeatedStratifiedKFold(n_splits=k, n_repeats=repeats, random_state=seed)
    return float(cross_val_score(make_pipe(), Z, s, cv=cv, scoring="roc_auc").mean())


rng = np.random.default_rng(config.SEED)
rows = []
for name, m in scopes.items():
    Z, s, c = Ltop[m], split[m], cls[m]
    k = min(5, int((s == 0).sum()), int((s == 1).sum()))
    obs = cv_auc(Z, s, k, 5, config.SEED)             # 11b ile ayni protokol (5 x k kat)
    null = np.empty(N_PERM)
    for b in range(N_PERM):
        sp = s.copy()
        for v in np.unique(c):                          # etiketleri SINIF ICINDE karistir
            idx = np.where(c == v)[0]
            sp[idx] = rng.permutation(s[idx])
        null[b] = cv_auc(Z, sp, k, 1, config.SEED + b)  # permutasyon basina 1 x k kat
    p = (1 + int((null >= obs).sum())) / (1 + N_PERM)
    rows.append({"kapsam": name, "gozlenen_AUC": obs, "bos_AUC_ort": float(null.mean()),
                 "bos_AUC_std": float(null.std()), "bos_AUC_%95": float(np.quantile(null, 0.95)),
                 "p_tek_tarafli": p, "n_perm": N_PERM})
    print(f"[{name}] gozlenen AUC={obs:.3f} | bos dagilim ort={null.mean():.3f} "
          f"std={null.std():.3f} %95={np.quantile(null, 0.95):.3f} | p={p:.4f} (n_perm={N_PERM})")

config.FIG_DIR.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(config.FIG_DIR / "eda_golub_split_auc_permutation.csv", index=False)
print("Tablo ->", config.FIG_DIR / "eda_golub_split_auc_permutation.csv")
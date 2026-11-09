"""EDA adim 3: teknik/batch fark kontrolleri (Golub ve Armstrong).

1) Anlamli genlerin yonu: AML'de yuksek / dusuk orani (dengeli mi, tek yonlu mu?)
2) Genel ifade duzeyi kaymasi: ornek basina medyan log2 (sinifa gore)
3) Kontrol probu testi: AFFX- (spike-in) problari sinifla degisiyor mu? (teknik fark gostergesi)
4) Golub: egitim/test farki, sinif icinde (sinif bilesimi etkisini ayirmak icin)
5) Armstrong: ALL vs MLL farki, AML vs (ALL+MLL) farkiyla karsilastirma

Yalnizca ogrenme amaclidir. Etiket kullanan testler modelleme/degerlendirmede kullanilmaz.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config, eda
from src.armstrong_loading import binary_labels, load_armstrong
from src.data_loading import load_golub

SPIKE = ["BioB", "BioC", "BioD", "Cre", "Dap", "Lys", "Phe", "Thr", "Trp"]


def summarize_de(name, L, y, ids):
    de = eda.diff_expression(L, y, ids, ids)
    sig = de[de["fdr"] < 0.05]
    up = (sig["log2fc"] > 0).mean() if len(sig) else float("nan")
    print(f"[{name}] test edilen {len(de)}, FDR<0.05: {len(sig)} (%{100 * len(sig) / max(len(de), 1):.1f}); "
          f"anlamli genlerin %{100 * up:.1f}'i AML'de YUKSEK, %{100 * (1 - up):.1f}'i AML'de DUSUK")
    return de


def median_shift(name, L, y):
    med = np.median(L, axis=1)
    a, b = med[y == 0], med[y == 1]
    p = mannwhitneyu(a, b).pvalue
    print(f"[{name}] ornek medyani (log2): ALL {a.mean():.3f}+-{a.std():.3f}, AML {b.mean():.3f}+-{b.std():.3f}, Mann-Whitney p = {p:.3g}")


# ---------------- Golub ----------------
X_tr, y_tr, X_te, y_te = load_golub()
raw_ids = pd.read_csv(config.TRAIN_FILE, usecols=["Gene Accession Number"])["Gene Accession Number"].astype(str).to_numpy()
X_all = pd.concat([X_tr, X_te]).to_numpy()
y_all = pd.concat([y_tr, y_te]).to_numpy()
split = np.array(["train"] * len(X_tr) + ["test"] * len(X_te))
is_ctrl = np.array([i.startswith("AFFX-") for i in raw_ids])
dm = eda.dudoit_mask(X_all[:, ~is_ctrl])
L = eda.to_log2(X_all[:, ~is_ctrl])[:, dm]
ids = raw_ids[~is_ctrl][dm]

print("===== GOLUB =====")
summarize_de("Golub", L, y_all, ids)
median_shift("Golub", L, y_all)
med = np.median(L, axis=1)
for cname, c in [("ALL", 0), ("AML", 1)]:
    m = y_all == c
    a, b = med[m & (split == "train")], med[m & (split == "test")]
    p = mannwhitneyu(a, b).pvalue
    print(f"[Golub, sinif ici] {cname}: egitim n={len(a)} vs test n={len(b)}, ornek medyani Mann-Whitney p = {p:.3g}")
ctrl_ids = raw_ids[is_ctrl]
Lc = eda.to_log2(X_all[:, is_ctrl])
spk = np.array([any(s in i for s in SPIKE) for i in ctrl_ids])
if spk.sum() > 3:
    de = eda.diff_expression(Lc[:, spk], y_all, ctrl_ids[spk], ctrl_ids[spk])
    print(f"[Golub] spike-in kontrol probu: {len(de)} test, FDR<0.05: {int((de['fdr'] < 0.05).sum())}")

# ---------------- Armstrong ----------------
XA, y3 = load_armstrong()
yb = binary_labels(y3, mll="ALL").to_numpy()
XC, _ = load_armstrong(drop_controls=False)
dmA = eda.dudoit_mask(XA.to_numpy())
LA = eda.to_log2(XA.to_numpy())[:, dmA]
idsA = np.array(XA.columns)[dmA]

print("\n===== ARMSTRONG =====")
summarize_de("Armstrong, AML vs (ALL+MLL)", LA, yb, idsA)
median_shift("Armstrong", LA, yb)
ctrl_cols = [c for c in XC.columns if c.startswith("AFFX-")]
spkA = [c for c in ctrl_cols if any(s in c for s in SPIKE)]
if len(spkA) > 3:
    LcA = eda.to_log2(XC[spkA].to_numpy())
    de = eda.diff_expression(LcA, yb, np.array(spkA), np.array(spkA))
    n_sig = int((de["fdr"] < 0.05).sum())
    print(f"[Armstrong] spike-in kontrol probu: {len(de)} test, FDR<0.05: {n_sig}")
    print(de.head(5)[["probe", "t", "fdr", "log2fc"]].round(4).to_string(index=False))

# ALL vs MLL (AML disarida), ayni gen kumesi
m = y3.to_numpy() != "AML"
yy = (y3.to_numpy()[m] == "MLL").astype(int)
de2 = eda.diff_expression(LA[m], yy, idsA, idsA)
print(f"[Armstrong] ALL vs MLL (n={int((yy == 0).sum())} vs {int((yy == 1).sum())}): "
      f"test {len(de2)}, FDR<0.05: {int((de2['fdr'] < 0.05).sum())}")
print("Not: AML vs ALL+MLL (n=28 vs 44) ile ALL vs MLL (n=24 vs 20) orneklem boyutlari farkli; "
      "karsilastirma yalnizca kabaca yol gostericidir.")

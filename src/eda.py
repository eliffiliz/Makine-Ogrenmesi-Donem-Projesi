"""EDA yardimci fonksiyonlari.

DIKKAT: Bu fonksiyonlar yalnizca KESIF amaclidir. Burada tum ornekler uzerinde yapilan
gen secimi / olcekleme, model egitimi ve degerlendirmesinde KULLANILMAZ (sizinti olur).
Modelleme asamasinda ayni adimlar yalnizca egitim katlarinda ogrenilecektir.
"""
import gzip
import io

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

from . import config

FLOOR, CEIL = 100.0, 16000.0  # sabit esikler (veriden ogrenilmez)


def to_log2(X):
    """Sabit esiklerle kirp, sonra log2 al."""
    return np.log2(np.clip(np.asarray(X, dtype=float), FLOOR, CEIL))


def dudoit_mask(X_raw, floor=100.0, ceil=16000.0):
    """Dudoit ve ark. (2002) gen filtresi (etiket kullanmaz).

    Kirpilmis [100, 16000] ham degerlerde max/min > 5 VE max-min > 500 olan genler tutulur.
    Hemen hemen sabit (cogu orneklerde taban degerinde) genleri eler.
    """
    Xc = np.clip(np.asarray(X_raw, dtype=float), floor, ceil)
    mx, mn = Xc.max(axis=0), Xc.min(axis=0)
    return (mx / mn > 5) & (mx - mn > 500)


def top_variance(L, k=2000):
    """En yuksek varyansli k gen sutununun indeksleri (etiket kullanmaz)."""
    v = L.var(axis=0)
    k = min(k, int((v > 0).sum()))
    return np.argsort(v)[::-1][:k]


def zscore(M):
    sd = M.std(axis=0, ddof=0)
    sd = np.where(sd == 0, 1.0, sd)
    return (M - M.mean(axis=0)) / sd


def run_pca(M, n=10):
    p = PCA(n_components=min(n, M.shape[0] - 1), random_state=config.SEED)
    scores = p.fit_transform(M)
    return scores, p.explained_variance_ratio_


def bh_fdr(p):
    """Benjamini-Hochberg FDR duzeltmesi."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def diff_expression(L, y, probe_ids, symbols):
    """Her gen icin Welch t-testi (AML=1 ile ALL=0), BH-FDR ve log2 kat degisimi.

    L: (ornek x gen) log2 matris, y: 0/1 etiket. Donus: p-degerine gore sirali DataFrame.
    """
    y = np.asarray(y)
    a, b = L[y == 0], L[y == 1]
    keep = L.std(axis=0) > 0
    t, p = stats.ttest_ind(b[:, keep], a[:, keep], equal_var=False)
    ok = ~np.isnan(p)
    idx = np.where(keep)[0][ok]
    t, p = t[ok], p[ok]
    df = pd.DataFrame({
        "probe": np.asarray(probe_ids)[idx],
        "symbol": np.asarray(symbols)[idx],
        "t": t,
        "p_value": p,
        "fdr": bh_fdr(p),
        "log2fc": b[:, idx].mean(axis=0) - a[:, idx].mean(axis=0),
        "col_index": idx,
    })
    return df.sort_values("p_value").reset_index(drop=True)


def probe_symbols(gpl, ids):
    """Onbellekteki GEO platform tablosundan prob -> ilk gen sembolu. Dosya yoksa prob adini dondurur."""
    path = config.DATA_DIR / "annotation" / f"{gpl}_family.soft.gz"
    if not path.exists():
        return list(ids)
    lines, inside = [], False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                inside = True
            elif line.startswith("!platform_table_end"):
                break
            elif inside:
                lines.append(line)
    tab = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
    col = next((c for c in ["Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol"] if c in tab.columns), None)
    if col is None:
        return list(ids)
    first = {}
    for pid, sym in zip(tab["ID"], tab[col]):
        if isinstance(sym, str):
            parts = [s.strip() for s in sym.split("///") if s.strip() and s.strip() != "---"]
            if parts:
                first[pid] = parts[0]
    return [first.get(i, i) for i in ids]

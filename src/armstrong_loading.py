"""Armstrong ve ark. (2002) MLL losemi verisini (MLL.arff) yukler.

ARFF dosyasinda ornek kimligi YOKTUR; satir sirasi sinifa gore dizilidir
(ALL, MLL, AML sirasiyla). Bu yuzden bolme/CV yaparken satir sirasi kullanilmaz,
her zaman karistirilarak ve sinifa gore tabakali bolunur.
"""
import numpy as np
import pandas as pd

from . import config

ARFF_PATH = config.DATA_DIR / "armstrong" / "raw" / "MLL.arff"
EXPECTED_COUNTS = {"ALL": 24, "MLL": 20, "AML": 28}


def _attr_name(line):
    return line.split(None, 2)[1].strip("'\"")


def load_armstrong(path=None, drop_controls=True, validate=True):
    """Donus: X (ornek x prob, float), y (ALL / MLL / AML metin etiketi).

    drop_controls=True: 'AFFX-' ile baslayan Affymetrix kontrol problarini atar.
    """
    path = path or ARFF_PATH
    names, rows, in_data = [], [], False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("%"):
                continue
            low = s.lower()
            if in_data:
                rows.append(s.split(","))
            elif low.startswith("@attribute"):
                names.append(_attr_name(s))
            elif low.startswith("@data"):
                in_data = True

    if names[-1].lower() != "class":
        raise ValueError(f"Son ozellik 'class' olmali, bulundu: {names[-1]}")
    probes = names[:-1]
    if len(set(probes)) != len(probes):
        raise ValueError("Tekrar eden prob adlari var")

    X = pd.DataFrame(
        np.array([[float(v) for v in r[:-1]] for r in rows]),
        columns=probes,
        index=[f"ARM_{i + 1:03d}" for i in range(len(rows))],  # yapay kimlik: sira numarasi
    )
    y = pd.Series([r[-1].strip() for r in rows], index=X.index, name="class")

    n_controls = int(sum(p.startswith("AFFX-") for p in probes))
    if drop_controls:
        X = X.loc[:, [not p.startswith("AFFX-") for p in X.columns]]

    if validate:
        counts = y.value_counts().to_dict()
        assert counts == EXPECTED_COUNTS, f"Sinif sayilari beklenenden farkli: {counts}"
        assert X.shape[0] == 72, f"Ornek sayisi 72 olmali: {X.shape[0]}"
        assert not np.isnan(X.to_numpy()).any(), "NaN deger var"
    X.attrs["n_control_probes"] = n_controls
    return X, y


def binary_labels(y, mll="ALL"):
    """ALL=0, AML=1. mll='ALL': MLL ornekleri ALL sayilir (n=72); mll='drop': MLL cikarilir (n=52)."""
    if mll == "ALL":
        return y.map({"ALL": 0, "MLL": 0, "AML": 1}).astype(int)
    if mll == "drop":
        return y[y != "MLL"].map({"ALL": 0, "AML": 1}).astype(int)
    raise ValueError("mll 'ALL' veya 'drop' olmali")

"""Golub ALL/AML verisini yükler ve etiketleri hasta numarası üzerinden eşleştirir.

Önemli: CSV'lerde hasta sütunları numara sırasında DEĞİLDİR. Etiketler sütun
konumuna göre değil, hasta numarasına göre eşleştirilir ve doğrulanır.
"""
import pandas as pd

from . import config

GENE_COLS = ["Gene Description", "Gene Accession Number"]


def load_expression(path):
    """Gen x örnek CSV'sini örnek x gen DataFrame'ine çevirir.

    'call' (P/A/M) sütunları atılır. İndeks hasta numarası (int) olur.
    """
    raw = pd.read_csv(path)
    # pandas tekrar eden 'call' başlıklarını 'call', 'call.1', ... yapar
    call_cols = [c for c in raw.columns if str(c).startswith("call")]
    raw = raw.drop(columns=call_cols)

    gene_ids = raw["Gene Accession Number"].astype(str)
    if gene_ids.duplicated().any():
        # tekrar eden accession'ları satır konumuyla benzersizleştir
        gene_ids = gene_ids + "_" + pd.Series(range(len(raw))).astype(str)

    patients = [c for c in raw.columns if c not in GENE_COLS]
    X = raw[patients].T
    X.columns = gene_ids.values
    X.index = X.index.astype(int)
    X.index.name = "patient"
    X = X.apply(pd.to_numeric, errors="coerce")
    return X.sort_index()


def load_labels(path):
    """actual.csv -> hasta numarasına göre indekslenmiş etiket Series'i (ALL=0, AML=1)."""
    lab = pd.read_csv(path)
    lab.columns = [c.strip().lower() for c in lab.columns]
    lab = lab.set_index("patient")["cancer"].astype(str).str.strip()
    lab.index = lab.index.astype(int)
    return lab.map(config.LABEL_MAP).astype(int), lab


def align(X, y_all):
    """X'in her satırına hasta numarasıyla etiket eşler; eksik eşleşmede hata verir."""
    missing = X.index.difference(y_all.index)
    if len(missing) > 0:
        raise ValueError(f"Etiketi bulunamayan hastalar: {list(missing)}")
    return y_all.loc[X.index]


def load_golub(train_path=None, test_path=None, label_path=None, validate=True):
    """Dönüş: X_train, y_train, X_test, y_test (örnek x gen, ALL=0 / AML=1)."""
    train_path = train_path or config.TRAIN_FILE
    test_path = test_path or config.TEST_FILE
    label_path = label_path or config.LABEL_FILE

    X_train = load_expression(train_path)
    X_test = load_expression(test_path)
    y_all, _ = load_labels(label_path)

    y_train = align(X_train, y_all)
    y_test = align(X_test, y_all)

    if validate:
        validate_golub(X_train, y_train, X_test, y_test)
    return X_train, y_train, X_test, y_test


def validate_golub(X_train, y_train, X_test, y_test):
    """Beklenen boyut, çakışma ve sınıf sayılarını kontrol eder."""
    assert X_train.shape[0] == 38, f"Eğitim örnek sayısı 38 olmalı: {X_train.shape[0]}"
    assert X_test.shape[0] == 34, f"Test örnek sayısı 34 olmalı: {X_test.shape[0]}"
    assert X_train.shape[1] == X_test.shape[1], "Gen sayıları eşleşmiyor"
    assert not set(X_train.index) & set(X_test.index), "Eğitim/test hasta çakışması var"
    assert list(X_train.columns) == list(X_test.columns), "Gen sırası eğitim ve testte farklı"

    inv = {v: k for k, v in config.LABEL_MAP.items()}
    tr = y_train.map(inv).value_counts().to_dict()
    te = y_test.map(inv).value_counts().to_dict()
    assert tr == config.EXPECTED_TRAIN_COUNTS, f"Eğitim sınıf sayıları beklenenden farklı: {tr}"
    assert te == config.EXPECTED_TEST_COUNTS, f"Test sınıf sayıları beklenenden farklı: {te}"

"""Sentetik veriyle yükleyicinin sütun sırasından bağımsız eşleştirmesini test eder."""
import numpy as np
import pandas as pd

from src.data_loading import load_expression, load_labels, align


def _make_csv(tmp_path, patients, name):
    # Hasta sütunları bilerek numara sırasında değil; her hastanın ardından 'call' sütunu
    rng = np.random.default_rng(0)
    data = {"Gene Description": ["g1", "g2", "g3"],
            "Gene Accession Number": ["A1", "A2", "A3"]}
    for i, p in enumerate(patients):
        data[str(p)] = rng.normal(size=3) + p
        data["call" if i == 0 else f"call.{i}"] = ["P", "A", "M"]
    path = tmp_path / name
    pd.DataFrame(data).to_csv(path, index=False)
    return path


def test_labels_follow_patient_number_not_column_order(tmp_path):
    patients = [5, 1, 3]  # karışık sıra
    path = _make_csv(tmp_path, patients, "expr.csv")
    lab = pd.DataFrame({"patient": [1, 3, 5], "cancer": ["ALL", "AML", "ALL"]})
    lab_path = tmp_path / "actual.csv"
    lab.to_csv(lab_path, index=False)

    X = load_expression(path)
    y_all, _ = load_labels(lab_path)
    y = align(X, y_all)

    assert list(X.index) == [1, 3, 5]
    assert list(y.values) == [0, 1, 0]  # 1->ALL, 3->AML, 5->ALL
    assert X.shape == (3, 3)  # call sütunları atıldı

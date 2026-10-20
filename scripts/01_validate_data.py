"""İlk veri kalitesi kontrolü: boyut, etiket eşleşmesi, INF/NaN/negatif sayıları."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.data_loading import load_golub

X_train, y_train, X_test, y_test = load_golub()  # içinde doğrulama var

print("Eğitim:", X_train.shape, "| Test:", X_test.shape)
print("Eğitim sınıf sayıları (0=ALL, 1=AML):", y_train.value_counts().to_dict())
print("Test sınıf sayıları  (0=ALL, 1=AML):", y_test.value_counts().to_dict())

for name, X in [("eğitim", X_train), ("test", X_test)]:
    v = X.to_numpy(dtype=float)
    print(f"[{name}] NaN: {np.isnan(v).sum()}, INF: {np.isinf(v).sum()}, "
          f"negatif: {(v < 0).sum()}, min: {np.nanmin(v):.1f}, max: {np.nanmax(v):.1f}")
print("Tüm kontroller geçti.")

"""Armstrong verisini yukler, temel kontrolleri yapar ve A3 tablosu icin sayilari yazdirir."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.armstrong_loading import binary_labels, load_armstrong

X, y = load_armstrong()
print(f"Kontrol probu (AFFX-) atildi: {X.attrs['n_control_probes']}")
print(f"Kalan matris: {X.shape[0]} ornek x {X.shape[1]} prob")
print("Sinif sayilari:", y.value_counts().to_dict())

yb = binary_labels(y, mll="ALL")
print(f"\nIkili (MLL -> ALL): n = {len(yb)}, ALL = {(yb == 0).sum()}, AML = {(yb == 1).sum()}, "
      f"ALL orani = {(yb == 0).mean():.3f}")
yd = binary_labels(y, mll="drop")
print(f"Ikili (MLL cikarildi): n = {len(yd)}, ALL = {(yd == 0).sum()}, AML = {(yd == 1).sum()}, "
      f"ALL orani = {(yd == 0).mean():.3f}")

v = X.to_numpy()
print(f"\nDeger araligi: min = {v.min():.0f}, max = {v.max():.0f}, medyan = {np.median(v):.0f}")
print(f"Negatif oran: {(v < 0).mean():.3f}, 20'nin alti: {(v < 20).mean():.3f}, 16000'in ustu: {(v > 16000).mean():.4f}")
print("\nSatir sirasi (ilk 3 / son 3 etiket):", y.iloc[:3].tolist(), y.iloc[-3:].tolist())
print("Ornek kimligi yapay (ARFF'de yok):", X.index[0], "...", X.index[-1])

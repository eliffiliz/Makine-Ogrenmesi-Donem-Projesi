"""Armstrong (MLL) verisi icin ilk inceleme: zip'i acar, dosyalari ve bicimini listeler.

Beklenen konum: data/armstrong/MLL.zip  (veya data/armstrong/ icine elle koyulan dosyalar)
Bu betik VERIYI YORUMLAMAZ; yalnizca dosya yapisini gosterir (gen kimlikleri var mi,
orneklerin sirasi/etiketleri nerede, boyutlar). Ciktiyi karar vermek icin kullan.
"""
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config

AR = config.DATA_DIR / "armstrong"
RAW = AR / "raw"
TEXT_EXT = {".txt", ".csv", ".tsv", ".data", ".dat", ".res", ".gct", ".names", ".cls", ""}

for z in AR.glob("*.zip"):
    RAW.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(z) as zf:
        zf.extractall(RAW)
        print(f"Acildi: {z.name} -> {RAW}")

root = RAW if RAW.exists() else AR
files = sorted(p for p in root.rglob("*") if p.is_file())
if not files:
    sys.exit(f"Dosya bulunamadi: {root}")

print(f"\n{len(files)} dosya:")
for p in files:
    print(f"  {p.relative_to(root)}  ({p.stat().st_size:,} bayt)")

for p in files:
    ext = p.suffix.lower()
    print(f"\n=== {p.relative_to(root)} ===")
    if ext == ".mat":
        from scipy.io import loadmat
        m = loadmat(p)
        for k, v in m.items():
            if not k.startswith("__"):
                print(f"  degisken {k!r}: {getattr(v, 'shape', type(v))}")
    elif ext in TEXT_EXT:
        with open(p, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
        print(f"  satir sayisi: {len(lines)}")
        for ln in lines[:3]:
            cells = ln.rstrip("\n").replace(",", "\t").split() if "\t" not in ln else ln.rstrip("\n").split("\t")
            print(f"  [{len(cells)} sutun] {ln[:200].rstrip()}")
    else:
        print("  (okunamayan bicim, ad ve boyut yukarida)")

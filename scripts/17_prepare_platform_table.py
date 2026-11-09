"""GEO platform tablosunu (kucuk dosya) onbellek bicimine cevirir (67 GB'lik *_family.soft.gz YERINE).

Kullanim:
    python scripts/17_prepare_platform_table.py <indirilen_dosya> GPL570

<indirilen_dosya> su biçimlerden biri olabilir (duz metin veya .gz):
  - GEO platform sayfasindan "Download full table..." ile inen tablo (ilk satir 'ID<TAB>...' basligi)
  - GPL570.annot.gz (GEO annot dosyasi; '!platform_table_begin' isaretleri icerir)
Cikti: data/annotation/<GPL>_family.soft.gz  -> scripts/15_mile_eda.py ve 16_common_genes_mile.py bu dosyayi
dogrudan kullanir.
"""
import gzip
import io
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config

if len(sys.argv) != 3:
    sys.exit(__doc__)
src, gpl = Path(sys.argv[1]), sys.argv[2]
opener = gzip.open if src.suffix == ".gz" else open
with opener(src, "rt", encoding="utf-8", errors="replace") as f:
    lines = f.read().splitlines()

begin = next((i for i, l in enumerate(lines) if l.startswith("!platform_table_begin")), None)
if begin is not None:
    start = begin + 1
else:
    start = next((i for i, l in enumerate(lines) if l.startswith("ID\t") or l.startswith('"ID"\t')), None)
    if start is None:
        sys.exit("Tablo basligi ('ID<TAB>...') bulunamadi; dosya beklenen bicimde degil.")
end = next((i for i in range(start, len(lines)) if lines[i].startswith("!platform_table_end")), len(lines))
body = [l for l in lines[start:end] if l.strip() and not l.startswith(("#", "^", "!"))]
tab = pd.read_csv(io.StringIO("\n".join(body)), sep="\t", dtype=str)
sym = next((c for c in ["Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol"] if c in tab.columns), None)
print(f"{gpl}: {len(tab)} satir; sutunlar (ilk 8): {list(tab.columns)[:8]}")
if sym is None:
    sys.exit("Gen sembolu sutunu bulunamadi; sutunlar: " + ", ".join(tab.columns))
print(f"Sembol sutunu: {sym!r}; sembolu olan satir: {int(tab[sym].notna().sum())}")

out = config.DATA_DIR / "annotation" / f"{gpl}_family.soft.gz"
out.parent.mkdir(parents=True, exist_ok=True)
with gzip.open(out, "wt", encoding="utf-8") as f:
    f.write("!platform_table_begin\n")
    f.write(tab.to_csv(sep="\t", index=False))
    f.write("!platform_table_end\n")
print("Yazildi ->", out)

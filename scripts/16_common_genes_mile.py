"""Golub (GPL80) ve MILE (GPL570) icin ortak gen sembolu sayisi (A3: ortak gen sayisi).

GEO platform aciklama tablolarini data/annotation/ altinda onbellege alir. Indirme basarisiz olursa dosyayi tarayicidan indirip data/annotation/GPL570_family.soft.gz adiyla koyun.
Ciktilar: ekrana sayilar, data/common_genes_golub_mile.txt
"""
import gzip
import io
import sys
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config

ANN = config.DATA_DIR / "annotation"
ANN.mkdir(parents=True, exist_ok=True)
SYMBOL_COLS = ["Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol", "gene_symbol"]


def soft_url(gpl):
    n = gpl[3:]
    stub = "GPL" + (n[:-3] if len(n) > 3 else "") + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{stub}/{gpl}/soft/{gpl}_family.soft.gz"


def platform_table(gpl):
    path = ANN / f"{gpl}_family.soft.gz"
    if not path.exists():
        print(f"Indiriliyor: {soft_url(gpl)}")
        with urllib.request.urlopen(soft_url(gpl)) as r, open(path, "wb") as f:
            f.write(r.read())
    lines, inside = [], False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                inside = True
            elif line.startswith("!platform_table_end"):
                break
            elif inside:
                lines.append(line)
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)


def probe_to_symbols(gpl):
    tab = platform_table(gpl)
    col = next((c for c in SYMBOL_COLS if c in tab.columns), None)
    if col is None:
        sys.exit(f"{gpl}: gen sembolu sutunu bulunamadi. Mevcut sutunlar: {list(tab.columns)}")
    print(f"{gpl}: {len(tab)} satir, sembol sutunu = {col!r}")
    out = {}
    for pid, sym in zip(tab["ID"], tab[col]):
        if isinstance(sym, str):
            parts = [s.strip() for s in sym.split("///") if s.strip() and s.strip() != "---"]
            if parts:
                out[pid] = parts
    return out


def summarize(name, probe_ids, mapping):
    ids = [p for p in probe_ids if not p.startswith("AFFX-")]
    mapped = [p for p in ids if p in mapping]
    symbols = {s for p in mapped for s in mapping[p]}
    print(f"[{name}] prob (kontrol disi): {len(ids)}, esleme bulunan: {len(mapped)} "
          f"(%{100 * len(mapped) / max(len(ids), 1):.1f}), benzersiz gen sembolu: {len(symbols)}")
    return symbols


golub_ids = pd.read_csv(config.TRAIN_FILE, usecols=["Gene Accession Number"])["Gene Accession Number"].astype(str).tolist()
mile_ids = (config.DATA_DIR / "mile" / "mile_probes.txt").read_text(encoding="utf-8").split("\n")
print(f"Golub prob: {len(golub_ids)}; MILE prob (NaN olanlar cikarilmis): {len(mile_ids)}\n")

g80 = probe_to_symbols("GPL80")
g570 = probe_to_symbols("GPL570")
print()
s_golub = summarize("Golub (GPL80)", golub_ids, g80)
s_mile = summarize("MILE (GPL570)", mile_ids, g570)
common = sorted(s_golub & s_mile)
print(f"\nOrtak gen sembolu (Golub kesisim MILE): {len(common)}")
print(f"Golub sembollerinin %{100 * len(common) / len(s_golub):.1f}'i MILE'da; MILE sembollerinin %{100 * len(common) / len(s_mile):.1f}'i Golub'da")
print("Ilk 10 ortak gen:", common[:10])
(config.DATA_DIR / "common_genes_golub_mile.txt").write_text("\n".join(common), encoding="utf-8")
print("Kaydedildi -> data/common_genes_golub_mile.txt")

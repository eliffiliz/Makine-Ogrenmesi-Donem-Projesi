"""Golub (GPL80, Hu6800) ve Armstrong (GPL91 / GPL8300, U95A / U95Av2) icin ortak gen sembolu sayisi.

GEO'dan platform aciklama tablolarini indirir (data/annotation/ altinda onbellege alir),
prob -> gen sembolu eslemesi yapar ve iki veri setinin ortak gen sayisini hesaplar.
Indirme basarisiz olursa: dosyayi tarayicidan indirip data/annotation/ icine koy
(ad: GPL80_family.soft.gz gibi). Sutun adi bulunamazsa mevcut sutunlari yazdirip durur.
"""
import gzip
import io
import sys
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config
from src.armstrong_loading import load_armstrong

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
    """Donus: {prob_id: [sembol, ...]} (bos / '---' atilir, '///' ile ayrilan coklu semboller ayrilir)."""
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


# --- Golub probleri (ham CSV'den, ek sonek olmadan) ---
golub_ids = pd.read_csv(config.TRAIN_FILE, usecols=["Gene Accession Number"])["Gene Accession Number"].astype(str).tolist()
print(f"Golub toplam prob: {len(golub_ids)}")
# --- Armstrong probleri ---
X, _ = load_armstrong(drop_controls=False)
arm_ids = list(X.columns)
print(f"Armstrong toplam prob: {len(arm_ids)}\n")

g80 = probe_to_symbols("GPL80")
a91 = probe_to_symbols("GPL91")
a8300 = probe_to_symbols("GPL8300")
print()

s_golub = summarize("Golub (GPL80)", golub_ids, g80)
s_arm91 = summarize("Armstrong (GPL91, U95A)", arm_ids, a91)
s_arm83 = summarize("Armstrong (GPL8300, U95Av2)", arm_ids, a8300)

common83 = sorted(s_golub & s_arm83)
common91 = sorted(s_golub & s_arm91)
print(f"\nOrtak gen sembolu (Golub kesisim Armstrong/GPL8300): {len(common83)}")
print(f"Ortak gen sembolu (Golub kesisim Armstrong/GPL91):   {len(common91)}")
print(f"Iki Armstrong aciklamasi arasi sembol farki: {len(s_arm91 ^ s_arm83)}")
print("Ilk 10 ortak gen:", common83[:10])

(config.DATA_DIR / "common_genes_golub_armstrong.txt").write_text("\n".join(common83), encoding="utf-8")
print("Kaydedildi -> data/common_genes_golub_armstrong.txt")

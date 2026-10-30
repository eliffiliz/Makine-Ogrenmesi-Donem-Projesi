"""MILE (GSE13159): ALL vs AML etiketlerini uretir ve doku kaynagi (BM/PB) ile karisiklik kontrolu yapar.

Girdi : data/mile/mile_metadata.csv  (03_mile_metadata.py uretir)
Cikti : data/mile/mile_labels_all_aml.csv  (yalnizca ALL ve AML ornekleri; ALL=0, AML=1)
"""
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import chi2_contingency

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import config

SRC = config.DATA_DIR / "mile" / "mile_metadata.csv"
OUT = config.DATA_DIR / "mile" / "mile_labels_all_aml.csv"


def lineage(cls):
    """Alt tip etiketinden soy (ALL / AML / other) cikarir."""
    if cls.startswith("AML"):
        return "AML"
    if "ALL" in cls:  # T-ALL, c-ALL/Pre-B-ALL, Pro-B-ALL, ALL with ..., mature B-ALL
        return "ALL"
    return "other"  # CLL, CML, MDS, Non-leukemia and healthy bone marrow


df = pd.read_csv(SRC)
df = df.rename(columns={"char_leukemia class": "subtype", "char_sample type": "sample_type"})
df["lineage"] = df["subtype"].map(lineage)

print("Tum 2096 ornek, soy dagilimi:")
print(df["lineage"].value_counts().to_string(), "\n")
print("'other' icindeki siniflar (dislanacak):")
print(df.loc[df["lineage"] == "other", "subtype"].value_counts().to_string(), "\n")

sub = df[df["lineage"].isin(["ALL", "AML"])].copy()
sub["label"] = sub["lineage"].map({"ALL": 0, "AML": 1})
n = len(sub)
print(f"ALL vs AML alt kumesi: n = {n}")
print(sub["lineage"].value_counts().to_string())
print(f"ALL orani: {(sub['label'] == 0).mean():.3f}  (Golub: {47/72:.3f})\n")

print("Soy x doku kaynagi (adet):")
ct = pd.crosstab(sub["lineage"], sub["sample_type"])
print(ct.to_string(), "\n")
print("Soy x doku kaynagi (satir yuzdesi):")
print((pd.crosstab(sub["lineage"], sub["sample_type"], normalize="index") * 100).round(1).to_string(), "\n")

chi2, p, dof, _ = chi2_contingency(ct)
print(f"Ki-kare bagimsizlik testi: chi2 = {chi2:.2f}, dof = {dof}, p = {p:.3g}")
print("Alt tip x doku kaynagi:")
print(pd.crosstab(sub["subtype"], sub["sample_type"]).to_string(), "\n")

sub[["geo_accession", "lineage", "label", "subtype", "sample_type"]].to_csv(OUT, index=False)
print("Kaydedildi ->", OUT)

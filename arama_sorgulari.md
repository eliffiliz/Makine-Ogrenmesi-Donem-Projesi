# Literatür arama sorguları

Sorgular 7 Kasım 2026'da çalıştırılmıştır. Scopus: "Advanced document search"; IEEE Xplore: "Advanced Search → Command Search" (yıl filtresi yok). Google Scholar'da ayrıca serbest anahtar kelime aramaları yapılmış, sorgu başına kayıt sayısı tutulmamıştır.

## Kayıt sayıları

| Veritabanı | S1 | S2 | S3 | S4 | S5 | Toplam |
|---|---|---|---|---|---|---|
| IEEE Xplore | 22 | 20 | 248 | 25 | 4 | 319 |
| Scopus | 13* | 98 | 861 | 10 | 13 | 995 |

\* Yalnızca 2021–2026. Genel toplam: 1.314 kayıt.

## Scopus

```
S1: TITLE-ABS-KEY(("Golub" OR "ALL/AML" OR "ALL vs AML") AND leukemia AND (microarray OR "gene expression") AND (classification OR "feature selection")) AND PUBYEAR > 2020 AND PUBYEAR < 2027
S2: TITLE-ABS-KEY(("selection bias" OR "data leakage" OR "nested cross-validation") AND (microarray OR "gene expression") AND classification) AND PUBYEAR > 1998 AND PUBYEAR < 2027
S3: TITLE-ABS-KEY(("deep learning" OR transformer OR TabPFN OR autoencoder) AND "gene expression" AND (leukemia OR cancer) AND classification) AND PUBYEAR > 1998 AND PUBYEAR < 2027
S4: TITLE-ABS-KEY(("MILE" OR "GSE13159" OR "Microarray Innovations in Leukemia") AND (leukemia OR AML OR ALL) AND classification) AND PUBYEAR > 1998 AND PUBYEAR < 2027
S5: TITLE-ABS-KEY(("GSE9476" OR "CuMiDa") AND leukemia AND classification) AND PUBYEAR > 1998 AND PUBYEAR < 2027
```

## IEEE Xplore

```
S1: ("All Metadata":"Golub" OR "All Metadata":"ALL/AML" OR "All Metadata":"ALL vs AML") AND "All Metadata":leukemia AND ("All Metadata":microarray OR "All Metadata":"gene expression") AND ("All Metadata":classification OR "All Metadata":"feature selection")
S2: ("All Metadata":"selection bias" OR "All Metadata":"data leakage" OR "All Metadata":"nested cross-validation") AND ("All Metadata":microarray OR "All Metadata":"gene expression") AND "All Metadata":classification
S3: ("All Metadata":"deep learning" OR "All Metadata":transformer OR "All Metadata":TabPFN OR "All Metadata":autoencoder) AND "All Metadata":"gene expression" AND ("All Metadata":leukemia OR "All Metadata":cancer) AND "All Metadata":classification
S4: ("All Metadata":"MILE" OR "All Metadata":"GSE13159" OR "All Metadata":"Microarray Innovations in Leukemia") AND ("All Metadata":leukemia OR "All Metadata":AML OR "All Metadata":ALL) AND "All Metadata":classification
S5: ("All Metadata":"GSE9476" OR "All Metadata":"CuMiDa") AND "All Metadata":leukemia AND "All Metadata":classification
```

## Dahil etme / dışlama

Dahil: (i) mikrodizi lösemi sınıflandırması veya uygulanabilir doğrulama/öznitelik seçimi metodolojisi, (ii) doğrulama protokolünün raporlanması, (iii) doğrulanabilir DOI. Dışlama: (i) kan yayması/mikroskopi görüntüsü, (ii) yalnızca RNA-seq (temel model kaynakları hariç), (iii) tam metne erişilemeyen kayıtlar, (iv) tekrarlar. Başlık/özet incelemesi, kartopu yöntemi ve Scholar aramalarıyla 28 kaynak dahil edilmiştir.

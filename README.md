# ALL/AML Lösemi Alt Tip Sınıflandırması (Gen İfadesi Verisi)

Dönem projesi: mikrodizi gen ifadesi verisinden ALL vs. AML alt tip sınıflandırması.
Çalışma bireyseldir. Problem "ALL vs. AML alt tip sınıflandırması" olarak tanımlanır
(her ikisi de lösemi alt tipidir).

- Ara rapor PDF: _eklenecek_
- Tanıtım videosu: _URL eklenecek_

## Kurulum

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Veri

```bash
python scripts/00_download_data.py   # Kaggle API token gerekir, bkz. data/README.md
python scripts/01_validate_data.py   # boyut ve etiket eşleşmesi kontrolü
pytest                               # birim testleri
```

## Klasör yapısı

| Klasör | İçerik |
|---|---|
| `src/` | yükleme, ön işleme ve model kodu |
| `scripts/` | çalıştırılabilir adımlar (sıra numaralı) |
| `tests/` | birim testleri |
| `figures/` | koddan üretilen şekiller |
| `report/` | rapor kaynakları (LaTeX) |
| `data/` | ham veri (depoya konmaz) |

## Tekrarlanabilirlik

- Sabit seed: `src/config.py` içinde `SEED = 42`
- Etiketler hasta numarası üzerinden eşleştirilir (CSV'de sütun sırası karışıktır)
- Rapordaki her şekil ve tablo bu depodaki bir script çıktısıyla eşleşecektir

## GenAI kullanımı

GenAI araçlarının kullanımı ayrı bir beyan dosyasında kayıt altındadır (`GENAI_LOG.md`).

## Veri kaynaklari

| Veri seti | Kaynak | Konum |
|---|---|---|
| Golub (1999), ana veri seti | Kaggle "Gene expression dataset (Golub)", doi:10.1126/science.286.5439.531 | `data/golub/` |
| Armstrong (2002), ikinci veri seti | https://csse.szu.edu.cn/staff/zhuzx/MLL.zip (GEO erisim numarasi yok), doi:10.1038/ng765 | `data/armstrong/MLL.zip` |
| Platform aciklamalari | GEO GPL80, GPL91, GPL8300 (09 indirir) | `data/annotation/` |
| MILE (GSE13159), yedek aday | GEO (03-05 script'leri) | `data/mile/` |

Ham veri depoya konmaz (`.gitignore`).

## EDA'yi yeniden uretme

```bash
python scripts/run_eda.py
```

| Script | Ne yapar |
|---|---|
| 01 | Golub'u yukler, etiketleri hasta numarasiyla eslestirir, boyut/sinif kontrolleri |
| 02 | Golub veri kalitesi, log oncesi/sonrasi dagilimlar |
| 03-05 | MILE metadata, ALL/AML etiketleri, ifade matrisi kontrolu (yedek aday) |
| 06-08 | Armstrong zip/ARFF inceleme, yukleyici ve dogrulama |
| 09 | Golub-Armstrong ortak gen sembolleri (GEO platform tablolari) |
| 10 | PCA, t-SNE, UMAP, ayirt edici genler, isi haritasi |
| 11 | Teknik/batch fark kontrolleri (spike-in, yon, egitim/test) |
| 12 | Armstrong veri kalitesi, birlesik Golub-Armstrong PCA |
| 13 | `call` sutunlari ve etiketsiz hiyerarsik kumeleme |

Sabit seed: `src/config.py` icinde `SEED = 42`. Kutuphane surumleri: `requirements-lock.txt`.

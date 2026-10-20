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

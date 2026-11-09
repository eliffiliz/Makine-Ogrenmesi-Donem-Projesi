# ALL/AML Lösemi Alt Tipi Sınıflandırması (Gen İfadesi Verisi)

Makine Öğrenmesi dönem projesi: mikrodizi gen ifadesi verisinden ALL vs. AML alt tip sınıflandırması.
Çalışma bireyseldir. Problem "ALL vs. AML alt tip sınıflandırması" olarak tanımlanır (her ikisi de lösemi alt tipidir;
modüldeki "non-cancer" ifadesi biyolojik olarak yanlıştır ve kullanılmaz).

| | |
|---|---|
| Öğrenci | [Elif Filiz · 503020260013] |
| Ara rapor (PDF) | `report/ara_rapor.pdf` (kaynak: `report/ara_rapor.tex`) |
| Tanıtım videosu | [video URL'si] |
| Aşama | Ara rapor: EDA ve veri analitiği (model eğitimi final aşamasındadır) |

## Araştırma soruları

1. **RQ1:** Golub'da 38/34 sabit bölme sonuçları ile 72 örnek üzerinde iç içe tabakalı çapraz doğrulama sonuçları arasında ne kadar fark vardır?
2. **RQ2:** Öznitelik seçiminin çapraz doğrulamadan önce yapılması (sızıntılı protokol), sızıntısız protokole göre doğruluğu ne kadar şişirir?
3. **RQ3:** Golub'da sızıntısız protokolle seçilen yöntemler, benzer bağımsız veri setlerinde (Armstrong ve MILE'ın ALL/AML alt kümesi) de benzer performans gösterir mi?

Ara raporda RQ1–RQ3 için gerekli EDA yapılmıştır (özellikle Golub eğitim/test ayrılabilirliği, ölçek farkları, ALL/AML ayrımının zorluğu). Modelleme sonuçları finaldedir.

## Kurulum

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt    # birebir sürümler için: pip install -r requirements-lock.txt
```

Araçlar: Python (numpy, pandas, scipy, scikit-learn, matplotlib; ayrıntı `requirements.txt`).

## Veri setleri

| Veri seti | Rol | Örnek × öznitelik | Sınıflar | Platform | Kaynak ve erişim |
|---|---|---|---|---|---|
| **Golub (1999)** | Ana veri seti | 72 × 7.129 | ALL 47 · AML 25 (eğitim 38, bağımsız test 34) | Affymetrix Hu6800 (GEO GPL80) | Kaggle "Gene expression dataset (Golub)"; doi:10.1126/science.286.5439.531. GEO kaydı bulunamadı. |
| **Armstrong (2002)** | Benzer veri seti | 72 × 12.533 | ALL 24 · MLL 20 · AML 28 | Affymetrix U95A / U95Av2 (GEO GPL91, GPL8300) | https://csse.szu.edu.cn/staff/zhuzx/MLL.zip; doi:10.1038/ng765. Özgün veri Broad Institute yayın sayfasında. GEO kaydı bulunamadı. |
| **MILE (GSE13159)** | Benzer veri seti | 1.292 × 54.630 (yalnızca ALL/AML) | ALL 750 · AML 542 | Affymetrix HG-U133 Plus 2.0 (GEO GPL570) | GEO **GSE13159** |

Golub ile ortak gen sembolü sayısı: Armstrong 5.822 · MILE 6.032 (`scripts/09_common_genes.py` ve `scripts/16_common_genes_mile.py` çıktısı).

### Lisans ve kullanım koşulları

| Veri seti | Lisans / koşul |
|---|---|
| Golub (Kaggle kopyası) | Kaggle veri seti sayfasında lisans **CC0: Public Domain** olarak görünür (teslimden önce canlı sayfadan yeniden kontrol edilmelidir). CC0 bildirimini Kaggle'daki yükleyici yapmıştır; özgün veri Golub ve ark. (1999)'a aittir, atıf yapılır. |
| Armstrong | Açık bir lisans veya kullanım koşulu bulunamadı. Özgün veri Broad Institute yayın sayfasında yer alır; indirilen `MLL.zip` bir üçüncü taraf (Shenzhen Üniversitesi) kopyasıdır. Veri yalnızca akademik/ders amacıyla kullanılır, Armstrong ve ark. (2002)'ye atıf yapılır, ham veri depoya konmaz. |
| MILE (GEO) | NCBI, GEO verisinin kullanımına veya dağıtımına kısıt koymaz; ancak bazı veri gönderenler telif/fikri mülkiyet hakkı iddia edebilir (GEO Disclaimer: https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html). GSE13159 ve Haferlach ve ark. (2010)'a atıf yapılır. |
| Bu depodaki kod | MIT (bkz. `LICENSE`). |

Her veri setinde özgün yayına atıf yapılır. Ham veri depoya konmaz (`.gitignore`).

## Veriyi indirme

```bash
python scripts/00_download_data.py   # Golub (Kaggle API token gerekir, bkz. data/README.md)
python scripts/01_validate_data.py   # boyut ve etiket eşleşmesi kontrolü
pytest                               # birim testleri
```

- **Armstrong:** `MLL.zip` yukarıdaki adresten elle indirilip `data/armstrong/MLL.zip` olarak konur (indirme scripti yok).
- **MILE:** `scripts/03_mile_metadata.py` GEO'dan `GSE13159_series_matrix.txt.gz` dosyasını indirir.
- **Platform tabloları:** `scripts/09_common_genes.py` GEO'dan GPL80, GPL91 ve GPL8300 tablolarını çeker (internet gerekir).

## EDA'yı yeniden üretme

```bash
python scripts/run_eda.py            # Golub + Armstrong EDA (01, 02, 06, 08, 09, 10, 11, 12, 13)
```

`run_eda.py` MILE ve eğitim/test ayırıcı analizlerini **çalıştırmaz**; bunlar ayrı çalıştırılır:

```bash
python scripts/03_mile_metadata.py && python scripts/04_mile_labels.py && python scripts/05_mile_expression.py
python scripts/15_mile_eda.py        # MILE EDA şekilleri
python scripts/16_common_genes_mile.py
python scripts/11c_split_auc_permutation.py   # eğitim–test ayırıcı AUC ve permütasyon testi
python scripts/14_split_distinguishability.py
```

| Script | Ne yapar |
|---|---|
| 01 | Golub'u yükler, etiketleri hasta numarasıyla eşleştirir, boyut ve sınıf kontrolleri |
| 02 | Golub veri kalitesi (eksik/INF, negatif ve uç değerler), log öncesi/sonrası dağılımlar |
| 03–05 | MILE metadata, ALL/AML etiketleri, ifade matrisi kontrolü (yedek aday) |
| 06–08 | Armstrong zip/ARFF inceleme, yükleyici ve doğrulama |
| 09 | Golub–Armstrong ortak gen sembolleri (GEO platform tabloları) |
| 10 | PCA (açıklanan varyans), t-SNE, UMAP, ayırt edici genler, ısı haritası |
| 11, 11b, 11c, 14 | Teknik/batch farkı kontrolleri (spike-in, yön, eğitim/test, ayırıcı AUC) |
| 12 | Armstrong veri kalitesi, birleşik Golub–Armstrong PCA |
| 13 | `call` (P/A/M) sütunları ve etiketsiz hiyerarşik kümeleme |
| 15–17 | MILE EDA, MILE ortak genler, platform tablosu hazırlama |

Şekiller ve tablolar `figures/` altına yazılır. Her şekil ve tablo bir script çıktısıyla eşleşir.

## Ana bulgular (EDA)

1. Golub'da etiketler hasta numarasıyla eşleştirilir; CSV'deki sütun sırası numara sırasında değildir. Sınıf oranı ≈ %65 ALL / %35 AML; ALL oranı eğitimde %71 (27/38), testte %59 (20/34).
2. Değerlerin yaklaşık %28–29'u negatif, %46–48'i 100'ün altında. 100'ün altındaki değerlerin ≈ %97'si `Absent` çağrılı. Ön işleme önerisi: [100; 16.000] aralığına kırp, log2 al.
3. PCA'da PC1 ve PC2 sırasıyla %13,2 ve %11,4 varyans açıklar; UMAP ve t-SNE'de ALL/AML daha net ayrışır.
4. Ayırt edici genler B hücre (CD19, CD79A/B, VPREB1) ve miyeloid (MPO, CD33, AZU1) belirteçleridir; CD79B ve TCL1A hem Golub hem Armstrong üst genleri arasındadır.
5. Eğitim ve test setini ayıran sınıflandırıcı AUC ≈ 0,78 alır (permütasyon p ≈ 0,005): hafif ama anlamlı bir batch kayması var (RQ1).
6. Golub ile Armstrong genel z-skorla birleştirildiğinde PC1 (%38,3) veri setini ayırır; veri seti başına z-skorla setler karışır ve sınıf yapısı korunur.
7. MILE'de PCA silhouette sınıfa göre 0,52, doku kaynağına göre −0,07; etiketsiz hiyerarşik kümelemede ARI 0,86.

Modelleme çıkarımları: sızıntısız pipeline (ölçekleme ve gen seçimi yalnızca CV katları içinde), stratified bölme ve sınıf ağırlığı, Balanced Accuracy/F1/ROC-AUC, 38/34 bölme yanında 72 örnekte nested CV.

## Tekrarlanabilirlik

- Sabit seed: `src/config.py` içinde `SEED = 42`.
- Kütüphane sürümleri: `requirements-lock.txt`.
- Etiketler hasta numarası üzerinden eşleştirilir (CSV'de sütun sırası karışıktır); `tests/` bunu doğrular.

## Klasör yapısı

| Klasör | İçerik |
|---|---|
| `src/` | yükleme, ön işleme ve model kodu |
| `scripts/` | çalıştırılabilir adımlar (sıra numaralı) |
| `tests/` | birim testleri |
| `figures/` | koddan üretilen şekiller ve özet tablolar |
| `report/` | ara rapor kaynağı (`ara_rapor.tex`; şekiller `../figures/` altından alınır) |
| `docs/` | literatür arama sorguları (`arama_sorgulari.md`) |
| `data/` | ham veri (depoya konmaz) |

## Finale kadar plan

- **Planlanan yöntemler:** L1 düzenlileştirmeli lojistik regresyon (gömülü gen seçimi, yorumlanabilirlik), SVM (küçük örnekli yüksek boyutlu veride güçlü referans), rastgele orman (doğrusal olmayan, ölçeklemeye duyarsız). Zaman kalırsa KNN ve XGBoost. Karar ağacı temel referans olarak raporlanır.
- **Doğrulama:** ön işleme, filtre, ölçekleme ve gen seçimi yalnızca eğitim katlarında; 5×5 tekrarlı iç içe tabakalı CV ve ek olarak 38/34 sabit bölme; en az 5 tohum; Accuracy, Balanced Accuracy, Precision, Recall, Specificity, F1, ROC-AUC, karışıklık matrisi; eşleştirilmiş istatistiksel test. Aynı protokol Armstrong ve MILE'a uygulanır.
- **Hedeflenen bonuslar:** B1 (MLP), B2 (TabPFN v2, öznitelik seçimi iç katlarda), B4 (Armstrong ve MILE birlikte); B3 yalnızca hocadan yazılı onay alınırsa.
- **Takvim:** 1. hf sızıntısız pipeline; 2. hf Golub'da üç yöntem ve iç içe CV; 3. hf MILE uyumlaştırması; 4. hf ablasyon ve yorumlanabilirlik; 5. hf bonuslar; 6. hf rapor ve video.
- **Riskler:** p ≫ n kaynaklı aşırı uyum ve yüksek varyans (tekrarlı CV, güven aralıkları); batch ve platform farkı; ortak gen eşlemesinde kayıp; MILE'da merkez etkisi; hesaplama süresi; Armstrong'un GEO numarasının olmaması.

## GenAI kullanımı

Kullanılan GenAI araçları ve doğrulama adımları `GENAI_LOG.md` dosyasında kayıt altındadır.

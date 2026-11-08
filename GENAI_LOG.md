### Yapay Zeka (GenAI) Kullanım Beyanı

Bu projede Yapay Zeka (Claude), kodlama taslakları oluşturma, hata ayıklama, görselleştirme düzeltmeleri, Grafik Yorumlama ve arama sorgusu formatlama aşamalarında yardımcı araç olarak kullanılmıştır. Projenin teorik çerçevesi, literatür taraması ve ana kararları tamamen araştırmacıya aittir.

| Süreç | GenAI Kullanıldı mı? | Kullanılan Araç | Geliştirici İpuçları (Prompts) | Değerlendirme & Uyarlamalar |
|---|---|---|---|---|
| **Literatür Taraması** | **Kısmen (Formatlama)** | Claude | *"Aşağıdaki anahtar kelimeleri akademik veritabanı arama formatına (IEEE/Scopus syntax) göre düzenle: ('Golub' OR 'ALL/AML' OR 'ALL vs AML') AND 'leukemia' AND ('microarray' OR 'gene expression') AND ('classification' OR 'feature selection')"* | Üretilen arama sorgusu akademik veritabanlarında literatür arama dizesi olarak kullanıldı. |
| **Fikir Geliştirme & Metodoloji** | **Kısmen** | Claude | *"Gen ifadesi verilerinde (p >> n) overfitting'i önlemek için hangi sınıflandırma ve regülasyon (L1/L2) stratejileri önerilir?"* | Metodolojik öneriler değerlendirildi, model seçimi ve cross-validation mimarisi tarafımdan kararlaştırıldı. |
| **Kodlama & Prototipleme** | **Evet** | Claude | *"Scikit-Learn ve Seaborn kullanarak PCA projeksiyonu ve Heatmap çizimi için temel fonksiyon şablonu (boilerplate code) öner."* | Üretilen kod taslağı incelendi, Golub ve Armstrong veri seti mimarisine göre parametreleri düzenlenip koda entegre edildi. |
| **Grafik & EDA Yorumlama** | **Evet** | Claude | *"PCA varyans açıklama oranı ve dendrogram çıktısını aldım. Bu grafik biyolojik/istatistiksel olarak ne anlama geliyor?"* | Claude'un Grafik ve EDA yorumları incelendi, doğruluğu literatürle teyit edilerek rapora aktarıldı. |
| **Görselleştirme & Hata Düzeltme** | **Evet** | Claude | *"Dağılım ve PCA grafiklerinde etiketler üst üste biniyor ve renk skalası sınıfları net ayırmıyor, Seaborn/Matplotlib kodunu nasıl düzeltirim?"* | Grafik eksen etiketlerinin okunaklı hale getirilmesi (rotasyon/layout düzenlemesi) ve renk paletlerinin düzeltilmesi sağlandı. |
| **Hata Ayıklama (Debugging)** | **Evet** | Claude | *"Özellik seçimi sonrası classifier'a veriyi sokarken 'ValueError: Shape mismatch' ve NaN hatası alıyorum, sebebi nedir?"* | Train/test split sonrası boyut uyumsuzlukları ve eksik veri imputasyon adımları Claude'a danışılarak giderildi. |
| **Metin Yazımı / Raporlama** | **Kısmen** | Claude | *"Hazırladığım teknik özetin dil bilgisi ve akademik akış kontrolünü yap."* | Çıktılar gözden geçirilerek düzenlendi ve rapora aktarıldı. |

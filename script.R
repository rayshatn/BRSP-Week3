#Analisis Ekspresi Gen Kanker Paru
#Dataset: GSE10072 (Lung Adenocarcinoma vs Normal)
#Platform: Microarray (Affymetrix GPL96)
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG) dan Enrichment 

#PERSIAPAN LINGKUNGAN KERJA (INSTALL & LOAD PACKAGE)
#1. Install BiocManager (manajer paket Bioconductor)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
#2. Install paket Bioconductor (GEOquery & limma)
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update =
                       FALSE)
#Install annotation package sesuai platform
#GPL96 = Affymetrix Human Genome U133A
BiocManager::install("hgu133a.db", ask = FALSE, update = FALSE)
#3. Install paket CRAN untuk visualisasi dan manipulasi data
install.packages(c("pheatmap", "ggplot2", "dplyr"))
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}
#4. Memanggil library
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133a.db)
library(AnnotationDbi)
library(umap)

#PENGAMBILAN DATA DARI GEO
#GSEMatrix = TRUE -> data diambil dalam format ExpressionSet
#AnnotGPL = TRUE -> anotasi gen (Gene Symbol) ikut diunduh
gset <- getGEO("GSE10072", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#PRE-PROCESSING DATA EKSPRESI
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm =
                            TRUE))
#LogTransform adalah variabel logika (TRUE / FALSE)
#Operator logika:
#> : lebih besar dari
#|| : OR (atau)
#&& : AND (dan)
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
#IF statement:
#Jika LogTransform = TRUE, maka lakukan log2
if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

#DEFINISI KELOMPOK SAMPEL
#pData(): metadata sampel
#source_name_ch1 berisi informasi kondisi biologis sampel
group_info <- pData(gset)[["source_name_ch1"]]
#make.names(): mengubah teks menjadi format valid untuk R
groups <- make.names(group_info)
#factor():mengubah data kategorik menjadi faktor
gset$group <- factor(groups)
#levels(): melihat kategori unik dalam faktor
nama_grup <- levels(gset$group)
print(nama_grup)

#DESIGN MATRIX (KERANGKA STATISTIK)
#model.matrix():membuat matriks desain untuk model linear
#~0 berarti TANPA intercept (best practice limma)
design <- model.matrix(~0 + gset$group)
#colnames(): memberi nama kolom agar mudah dibaca
colnames(design) <- levels(gset$group)
#Menentukan perbandingan biologis
grup_kanker <- nama_grup[1]
grup_normal <- nama_grup[2]
contrast_formula <- paste(grup_kanker, "-", grup_normal)
print(paste("Kontras yang dianalisis:", contrast_formula))

#ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)
#lmFit(): membangun model linear untuk setiap gen
fit <- lmFit(ex, design)
#makeContrasts(): mendefinisikan perbandingan antar grup
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels
                                 = design)
#contrasts.fit(): menerapkan kontras ke model
fit2 <- contrasts.fit(fit, contrast_matrix)
#eBayes():empirical Bayes untuk menstabilkan estimasi varians
fit2 <- eBayes(fit2)
#topTable():mengambil hasil akhir DEG
#adjust = "fdr" -> koreksi multiple testing
#p.value = 0.01 -> gen sangat signifikan
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)
head(topTableResults)

#ANOTASI NAMA GEN
#Mengambil ID probe dari hasil DEG
probe_ids <- rownames(topTableResults)
#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)
#Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)
topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)
#Cek hasil anotasi
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#BOXPLOT DISTRIBUSI NILAI EKSPRESI
group_colors <- as.numeric(gset$group)
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)
legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

#DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT)
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)
ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )

#UMAP (VISUALISASI DIMENSI RENDAH)
umap_input <- t(ex)
#Jalankan UMAP
umap_result <- umap(umap_input)
#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)
#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )

#VISUALISASI VOLCANO PLOT
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val <
                      0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val <
                      0.01] <- "DOWN"
#Visualisasi
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color =
                           status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" =
                                  "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Kanker Paru")

#VISUALISASI HEATMAP
#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]
top50 <- head(topTableResults, 50)
#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID, ]
#Gunakan Gene Symbol (fallback ke Probe ID)
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID, # jika SYMBOL kosong → probe ID
  top50$SYMBOL # jika ada → gene symbol
)
rownames(mat_heatmap) <- gene_label
#Pembersihan data (WAJIB agar tidak error hclust)
#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]
#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]
#Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)
#Visualisasi heatmap
pheatmap(
  mat_heatmap,
  scale = "row", # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE, # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

#ANALISIS GENE ONTOLOGY (GO) & KEGG
# Install package enrichment
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"), 
                     ask = FALSE, update = FALSE)
# Load library
library(clusterProfiler)
library(org.Hs.eg.db)
# Ambil gene symbol unik dan tidak NA
deg_genes <- unique(topTableResults$SYMBOL)
deg_genes <- deg_genes[!is.na(deg_genes)]

length(deg_genes)  # cek jumlah gen
#Konversi symbol menjadi entrez id
gene_df <- bitr(
  deg_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
entrez_ids <- gene_df$ENTREZID
#Analisis GO
ego <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",   # BP, MF, CC
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
head(ego)
#Visualisasi barplot dan dotplot
barplot(ego, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scale = "free")
dotplot(ego, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scale = "free")
#Analisis KEGG
ekegg <- enrichKEGG(
  gene         = entrez_ids,
  organism     = "hsa",   # Homo sapiens
  pvalueCutoff = 0.05
)
# Konversi agar gene symbol terbaca
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(ekegg)
#Visualisasi dotplot & barplot
dotplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment")
barplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment")

#MENYIMPAN HASIL
# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_GSE10072_DEG.csv")
write.csv(as.data.frame(ego), "GO_Enrichment_GSE10072.csv")
write.csv(as.data.frame(ekegg), "KEGG_Enrichment_GSE10072.csv")
message("Analisis selesai. File hasil telah disimpan.")
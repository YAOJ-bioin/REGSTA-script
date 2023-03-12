##########################################################################
rm(list=ls())
##########################################################################
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("Seurat")
#install.packages("SeuratObject")
#install.packages("dplyr")
#install.packages("ggthemes")
#install.packages("Matrix")
##########################################################################
##########################################################################
library(ggplot2)
library(ggpubr)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggthemes)
##########################################################################
##
##########################################################################
matrix_path <- "./data4test/matrix4test/"                
png_path <- "./data4test/scSTplot/"                                       
MT_path <- "./data4test/MT_CH_gene/rice_mt_genes_mini.txt"
CH_path <- "./data4test/MT_CH_gene/rice_pt_genes_mini.txt"
MiniCsv <- "./data4test/img4test_cp_outlines_scResGem_mini.csv"
ncbi_project_name <- "rice_embryo"
project_name <- "test"
##########################################################################
##
##########################################################################
raw_data <- Read10X(data.dir = matrix_path)  
obj <- CreateSeuratObject(counts = raw_data, project = "test", min.cells =3, min.features = 100)
##########################################################################
##QC
##########################################################################
MT_gene <- read.table(MT_path, header = F, sep = "\t")
MT_gene$V1 <- gsub("_","-", MT_gene$V1)
MT_gene <- MT_gene[MT_gene$V1 %in% rownames(obj),]
CH_gene <- read.table(CH_path, header = F, sep = "\t")
CH_gene$V1 <- gsub("_","-", CH_gene$V1)
CH_gene <- CH_gene[CH_gene$V1 %in% rownames(obj),]
obj[["percent.mt"]] <- PercentageFeatureSet(obj, features=MT_gene)
obj[["percent.ch"]] <- PercentageFeatureSet(obj, features = CH_gene)
QC01_path <- paste0(png_path, ncbi_project_name, "_", project_name, "_01beforeQCmetrics_100.png")
png(QC01_path, 1200, 900, res=150)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ch"), ncol = 4) & theme(axis.title = element_text(size = 11), axis.text = element_text(size = 8), plot.title = element_text(size=12)))
dev.off()
QC02_path <- paste0(png_path, ncbi_project_name, "_", project_name, "_02beforeQCrelation_100.png")
png(QC02_path, 2400, 900, res=150)
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend() + theme(axis.title = element_text(size = 11), axis.text = element_text(size = 8), plot.title = element_text(size=12))
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.ch") + NoLegend() + theme(axis.title = element_text(size = 11), axis.text = element_text(size = 8), plot.title = element_text(size=12))
plot3 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(axis.title = element_text(size = 11), axis.text = element_text(size = 8), plot.title = element_text(size=12))
print(plot1 + plot2 + plot3)
dev.off()
##########################################################################
##MT CH
##########################################################################
MT_gene <- read.table(MT_path, header = F, sep = "\t")
MT_gene <- MT_gene[MT_gene$V1 %in% rownames(obj),]
CH_gene <- read.table(CH_path, header = F, sep = "\t")
CH_gene <- CH_gene[CH_gene$V1 %in% rownames(obj),]
remove_gene <- c(MT_gene, CH_gene)
obj <- obj[!(rownames(obj) %in% remove_gene),]
##########################################################################
##########################################################################
MT_gene <- read.table(MT_path, header = F, sep = "\t")
MT_gene$V1 <- gsub("_","-", MT_gene$V1)
MT_gene <- MT_gene[MT_gene$V1 %in% rownames(obj),]
CH_gene <- read.table(CH_path, header = F, sep = "\t")
CH_gene$V1 <- gsub("_","-", CH_gene$V1)
CH_gene <- CH_gene[CH_gene$V1 %in% rownames(obj),]
obj[["percent.mt"]] <- PercentageFeatureSet(obj, features=MT_gene)
obj[["percent.ch"]] <- PercentageFeatureSet(obj, features = CH_gene)
obj <- subset(obj, subset = nFeature_RNA >= 100 & nFeature_RNA <= 5000 & percent.mt < 5)

QC03_path <- paste0(png_path, ncbi_project_name, "_", project_name, "_03_afterQCmetrics_100.png")
png(QC03_path, 1200, 900, res=150)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ch"), ncol = 4) & theme(axis.title = element_text(size = 11), axis.text = element_text(size = 8), plot.title = element_text(size=12)))
dev.off()
QC04_path <- paste0(png_path, ncbi_project_name, "_", project_name, "_04afterQCrelation_100.png")
png(QC04_path, 2400, 900, res=150)
plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend() + theme(axis.title = element_text(size = 11), axis.text = element_text(size = 8), plot.title = element_text(size=12))
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.ch") + NoLegend() + theme(axis.title = element_text(size = 11), axis.text = element_text(size = 8), plot.title = element_text(size=12))
plot3 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(axis.title = element_text(size = 11), axis.text = element_text(size = 8), plot.title = element_text(size=12))
print(plot1 + plot2 + plot3)
dev.off()
##########################################################################
##########################################################################
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- SCTransform(obj,  verbose = TRUE)
obj <- RunPCA(obj)
QC06_path <- paste0(png_path, ncbi_project_name, "_", project_name, "_06dimensionPick_100.png")
png(QC06_path, 2400, 900, res=150)
p1 <- VizDimLoadings(obj, dims = 1:2, reduction = "pca")
p2 <- DimPlot(obj, reduction = "pca")
ggarrange(p1, p2, ncol = 2, nrow = 1) 

QC07_path <- paste0(png_path, ncbi_project_name, "_", project_name, "_07dimensionPick_2_100.png")
png(QC07_path, 900, 900, res=150)
ElbowPlot(obj)
dev.off()

##########################################################################
##res=0.5
##########################################################################
### clustering
sce_0.5 <- FindNeighbors(obj, dims = 1:13)
sce_0.5 <- FindClusters(sce_0.5, resolution = 0.5)
sce_0.5 <- RunUMAP(sce_0.5, dims = 1:13)
sce_0.5 <- RunTSNE(sce_0.5, dims = 1:13 ,check_duplicates = FALSE)

QC08_path <- paste0(png_path, ncbi_project_name, "_", project_name, "_08umap_0.5_100.png")
png(QC08_path, 2400, 900, res=150)
p1 <- DimPlot(sce_0.5, 
              pt.size = 1,
              repel = T,
              label.size = 5,
              reduction = "umap",
              cols =  c("#DCB717","#E2272E","#5EA5C9","#76AB62","#91217F","#fdb462","#EA9014","#5A7FB9","#757C98","#178F3B",
                                 "#F1A1BF","#33B9C1","#B39BC2","#851B3B","#fb8072","#80b1d3","#b3de69","#fccde5","#d9d9d9","#bc80bd",
                                 "#9ebcda","#5d65d9"))
# label
p2 <- DimPlot(sce_0.5,
              pt.size = 1,
              repel = T,
              label.size = 5,
              reduction = "umap", label = TRUE,
              cols =  c("#DCB717","#E2272E","#5EA5C9","#76AB62","#91217F","#fdb462","#EA9014","#5A7FB9","#757C98","#178F3B",
                                 "#F1A1BF","#33B9C1","#B39BC2","#851B3B","#fb8072","#80b1d3","#b3de69","#fccde5","#d9d9d9","#bc80bd",
                                 "#9ebcda","#5d65d9"))
ggarrange(p1, p2, ncol = 2, nrow = 1) 
dev.off()

QC09_path <- paste0(png_path, ncbi_project_name, "_", project_name, "_09tsne_0.5_100.png")
png(QC09_path, 2400, 900, res=150)
p1 <- DimPlot(sce_0.5, 
              pt.size = 1,
              repel = T,
              label.size = 5,
              reduction = "tsne",
              cols =  c("#DCB717","#E2272E","#5EA5C9","#76AB62","#91217F","#fdb462","#EA9014","#5A7FB9","#757C98","#178F3B",
                                 "#F1A1BF","#33B9C1","#B39BC2","#851B3B","#fb8072","#80b1d3","#b3de69","#fccde5","#d9d9d9","#bc80bd",
                                 "#9ebcda","#5d65d9"))
# label
p2 <- DimPlot(sce_0.5, 
              pt.size = 1,
              repel = T,
              label.size = 5,
              reduction = "tsne", label = TRUE,
              cols =  c("#DCB717","#E2272E","#5EA5C9","#76AB62","#91217F","#fdb462","#EA9014","#5A7FB9","#757C98","#178F3B",
                                 "#F1A1BF","#33B9C1","#B39BC2","#851B3B","#fb8072","#80b1d3","#b3de69","#fccde5","#d9d9d9","#bc80bd",
                                 "#9ebcda","#5d65d9"))
ggarrange(p1, p2, ncol = 2, nrow = 1) 
dev.off()

sce_0.5.markers <- FindAllMarkers(sce_0.5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top6 <- sce_0.5.markers %>%  group_by(cluster) %>% slice_max(n = 6, order_by = avg_log2FC)

markers <- top6$gene
markers <- as.data.frame(markers)
markerdata <- ScaleData(sce_0.5, features = as.character(unique(markers$markers)), assay = "RNA")

QC15_path <- paste0(png_path, ncbi_project_name, "_", project_name, "_15heatmap_0.5_100.png")
png(QC15_path, 2400, 1800, res=150)
DoHeatmap(markerdata,
          features = as.character(unique(markers$markers)),
          group.by = "seurat_clusters",
          assay = 'RNA',
          group.colors = c("#91217F","#fdb462","#EA9014","#5A7FB9","#757C98","#178F3B","#F1A1BF","#33B9C1","#B39BC2","#851B3B","#fb8072","#80b1d3","#b3de69","#fccde5","#d9d9d9","#bc80bd",
                          "#9ebcda","#5d65d9"))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
dev.off()


#scST
miniCsv<- read.csv(MiniCsv,header = TRUE)
metadata<-sce_0.5@meta.data #
metadata$cell_id <-row.names(metadata) 
ST <- miniCsv %>%  
  mutate(cell_id = as.numeric(cell_id))
metadata <- metadata %>%  
  mutate(cell_id = as.numeric(cell_id))
#
ST <- ST %>% 
  left_join(metadata,by="cell_id")
ST <- ST[complete.cases(ST[, c(12)]), ]


QC16_path <- paste0(png_path, ncbi_project_name, "_", project_name, "_16scST_0.5_100.png")
png(QC16_path,3600, 4800, res=300)
ggplot(ST,mapping=aes(x=img_x ,y=img_y,color=seurat_clusters))+
  geom_point(size=0.01)+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(values = c("#76AB62","#178F3B","#5EA5C9","#91217F","#757C98","#EA9014","#fdb462","#E2272E","#851B3B","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9","#d9d9d9"),
                     #labels = c("RA","LS","EP-CR","SCL2","PL","COL2","COL2","SCL2")
                     )+
  theme_classic()+
  scale_y_continuous(breaks = NULL)+
  theme(panel.grid =element_blank()) +   
  theme(axis.text = element_blank()) +   
  theme(axis.ticks = element_blank()) +   
  theme(panel.border = element_blank())
dev.off()


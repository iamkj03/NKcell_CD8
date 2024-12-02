#Figure 3

library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(igraph)
library(readr)
library(data.table)
library(cowplot)
library(readxl)
table(cd8.combined$orig.ident)

data_dir_iso <- c("./CD8Iso/outs/")
data_dir_nk <- c("./CD8NK/outs/")
cd8NK <- Read10X( paste0(data_dir_nk,"filtered_feature_bc_matrix/") )
cd8iso <- Read10X( paste0(data_dir_iso,"filtered_feature_bc_matrix/") )

save(cd8NK,cd8iso, file="./beforeseurat_cd8.RData")
load("./beforeseurat_cd8.RData")

cd8NK <- CreateSeuratObject(counts = cd8NK, 
                            project = "cd8nk", 
                            min.cells = 3, #low quality genes
                            min.features = 200) #at least 200 features

cd8iso <- CreateSeuratObject(counts = cd8iso, 
                             project = "cd8iso", 
                             min.cells = 3, #low quality genes
                             min.features = 200) #at least 200 features



#qc-check and selecting cells -filter cells that have unique feature counts over 2,500 or less than 200, We filter cells that have >5% mitochondrial counts
cd8NK[["percent.mt"]] <- PercentageFeatureSet(cd8NK, pattern = "^mt-")
cd8iso[["percent.mt"]] <- PercentageFeatureSet(cd8iso, pattern = "^mt-")

VlnPlot(cd8NK, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(cd8iso, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1_nk <- FeatureScatter(cd8NK, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_nk <- FeatureScatter(cd8NK, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_nk + plot2_nk

plot1_iso <- FeatureScatter(cd8iso, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_iso <- FeatureScatter(cd8iso, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1_iso + plot2_iso

cd8iso <- subset(cd8iso, subset = nFeature_RNA < 7500)#nFeature_RNA > 1300)# & nFeature_RNA < 7500 )#& percent.mt < 15)
cd8NK <- subset(cd8NK, subset = nFeature_RNA < 7500 )#nFeature_RNA > 1400)# & nFeature_RNA < 7500 )#& percent.mt < 15)

#integrating
cd8 <- c(cd8iso,cd8NK)

cd8 <- lapply(X = cd8, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures((object.list = cd8))

cd8.anchors <- FindIntegrationAnchors(object.list = cd8, anchor.features = features)

# this command creates an 'integrated' data assay
cd8.combined <- IntegrateData(anchorset = cd8.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(cd8.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
cd8.combined <- ScaleData(cd8.combined, verbose = FALSE)
cd8.combined <- RunPCA(cd8.combined, npcs = 30, verbose = FALSE)
cd8.combined <- RunUMAP(cd8.combined, reduction = "pca", dims = 1:30)
cd8.combined <- FindNeighbors(cd8.combined, reduction = "pca", dims = 1:30)
cd8.combined <- FindClusters(cd8.combined, resolution = 0.5)

#save(cd8.combined,file="afternorm15_cd8inte_upperbound.RData")

#load("./afternorm_cd8inte_upperbound.RData")
load("afternorm15_cd8inte_upperbound.RData")

counts <- cd8.combined[["RNA"]]@counts
metadata <- cd8.combined@meta.data

save(counts,metadata,file = "./cd8_counts15_upperbound.RData")

#cell annotation
library(SingleR)
library(SingleCellExperiment)
library(scater)
library(tidyverse)

#load counts & metadata
#load('afternorm_cd8inte.RData')
sce <- SingleCellExperiment(assays = list(counts = counts))

# lognorm transform
sce <- logNormCounts(sce)

# load various reference databases
ref_imm <- ImmGenData()
ref_mouse <- MouseRNAseqData()


# cell type prediction
pred_imm <- SingleR(test = sce, ref = ref_imm, labels = ref_imm$label.main)
pred_mouse <- SingleR(test = sce, ref = ref_mouse, labels = ref_mouse$label.main)

pred_list <- list(pred_imm,
                  pred_mouse)

#upper
#save(pred_list,
#     file = "./cd8_norm_singleR_upperbound.RData")
save(pred_list,
     file = "./cd8_norm15_singleR_upperbound.RData")




#load("./cd8_norm_singleR_upperbound.RData")

#prediction diagnostics - HumanPrimaryCellAtlasData
pred_plot_imm <- plotScoreHeatmap(pred_list[[1]],
                                  annotation_col=as.data.frame(metadata[,"seurat_clusters",
                                                                        drop=FALSE]))
ggsave(plot = pred_plot_imm, 
       filename = "./cd8_pred_heatmap_imm15_upperbound.pdf", 
       device = "pdf",
       height  = 8,
       width  = 9)

pred_table_imm <- table(pred_list[[1]]$labels, metadata[,"seurat_clusters"])
write.csv(pred_table_imm, 
          file = "./cd8_pred_table_imm15_upperbound.csv",
          quote = FALSE,
          row.names = TRUE)


# prediction diagnostics - DatabaseImmuneCellExpressionData
pred_plot_mouse <- plotScoreHeatmap(pred_list[[2]], 
                                    annotation_col=as.data.frame(metadata[,"seurat_clusters",
                                                                          drop=FALSE]))
ggsave(plot = pred_plot_mouse, 
       filename = "./cd8_pred_heatmap_mouse15_upperbound.pdf", 
       device = "pdf",
       height  = 8,
       width  = 9)

pred_table_mouse <- table(pred_list[[2]]$labels, metadata[,"seurat_clusters"])
write.csv(pred_table_mouse, 
          file = "./cd8_pred_table_mouse15_upperbound.csv",
          quote = FALSE,
          row.names = TRUE)


# Visualization

pdf("./cd8.pdf")
p1 <- DimPlot(cd8.combined, reduction = "umap", group.by = "seurat_clusters",label=TRUE)
p2 <- DimPlot(cd8.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1
p2

dev.off()

Idents(cd8.combined) <- "orig.ident"
cd8.combined.nk <- subset(cd8.combined, idents = "cd8nk", invert = FALSE)
cd8.combined.iso <- subset(cd8.combined, idents = "cd8iso", invert = FALSE)

cd8.combined.nk
cd8.combined$nkiso_clu <- paste0(cd8.combined$orig.ident,"_",cd8.combined$seurat_clusters)

Idents(cd8.combined) <- "nkiso_clu"
cd8.combined.iso_1vsnk_5 <- subset(cd8.combined, idents=c("cd8iso_1","cd8nk_5"), invert=FALSE)
write.csv(cd8.combined.iso_1vsnk_5@assays$RNA@data,"./cd8_genematrix_iso_1vsnk_5.csv")
write.csv(cd8.combined.iso_1vsnk_5$seurat_clusters, "./cd8_label_iso_1vsnk_5.csv")

#two condition side by side
Idents(cd8.combined) <- "seurat_clusters"
Idents(cd8.combined.nk) <- "seurat_clusters"
Idents(cd8.combined.iso) <- "seurat_clusters"

cd8.combined <- subset(cd8.combined, idents = c(10,11), invert = TRUE)
cd8.combined.nk <- subset(cd8.combined.nk, idents = c(10,11), invert = TRUE)
cd8.combined.iso <- subset(cd8.combined.iso, idents = c(10,11), invert = TRUE)

cd8.combined.nk_1vs5 <- subset(cd8.combined.nk, idents=c(1,5), invert=FALSE)

write.csv(cd8.combined.nk_1vs5@assays$RNA@data,"./cd8nk_genematrix_1vs5.csv")
write.csv(cd8.combined.nk_1vs5$seurat_clusters, "./cd8nk_label_1vs5.csv")


p1 <- DimPlot(cd8.combined, reduction = "umap", split.by = "orig.ident",label=TRUE)
p2 <- DimPlot(cd8.combined.nk, reduction = "umap", split.by = "orig.ident",label=TRUE)
p3 <- DimPlot(cd8.combined.iso, reduction = "umap", split.by = "orig.ident",label=TRUE)
p1

abc <- which(cd8.combined.nk@reductions[["umap"]]@cell.embeddings[,1] >-8)
cd8.combined.nk <- cd8.combined.nk[,abc]

abc <- which(cd8.combined.iso@reductions[["umap"]]@cell.embeddings[,1] >-8)
cd8.combined.iso <- cd8.combined.iso[,abc]

abc <- which(cd8.combined@reductions[["umap"]]@cell.embeddings[,1] >-8)
cd8.combined <- cd8.combined[,abc]

p1 <- DimPlot(cd8.combined, reduction = "umap", split.by = "orig.ident",label=TRUE)
p2 <- DimPlot(cd8.combined.nk, reduction = "umap", split.by = "orig.ident",label=TRUE)
p3 <- DimPlot(cd8.combined.iso, reduction = "umap", split.by = "orig.ident",label=TRUE)


write.csv(cd8.combined.nk@assays$RNA@data,"./cd8nk_genematrix.csv")
write.csv(cd8.combined.nk$seurat_clusters, "./cd8nk_label.csv")

write.csv(cd8.combined.iso@assays$RNA@data,"./cd8iso_genematrix.csv")
write.csv(cd8.combined.iso$seurat_clusters, "./cd8iso_label.csv")

# Figure 3A
pdf("./Dimplot.pdf")

p2
p3
dev.off()

# Figure 3B
p1[[1]]
g <- ggplot_build(p1)
g
color_list <- as.vector(unique(g$data[[1]]["colour"]))

p1$data
Idents(cd8.combined) <- "seurat_clusters"

preprop_table <- table(cd8.combined$seurat_clusters,cd8.combined$orig.ident)
preprop_table
preprop_table[,1] <- preprop_table[,1]/sum(preprop_table[,1])
preprop_table[,2] <- preprop_table[,2]/sum(preprop_table[,2])
preprop_table

iso_pie <- preprop_table[1:10,1]
iso_pie <- as.data.frame(iso_pie)
iso_pie$group <- rownames(iso_pie)
nk_pie <- preprop_table[1:10,2]
nk_pie <- as.data.frame(nk_pie)
nk_pie$group <- rownames(nk_pie)

iso_pie$density <- 0.3
iso_pie$density[6] <- 1
iso_pie$density[2] <- 1

nk_pie$density <- 0.3
nk_pie$density[6] <- 1
nk_pie$density[2] <- 1

pdf("./pie_chart_new.pdf")
ggplot(iso_pie, aes(x = "", y = iso_pie, fill = group, alpha=density)) + scale_alpha_continuous(range=c(0,1),limits=c(0,1))+
  geom_col() +
  coord_polar(theta = "y")+ theme_void()+  guides(fill = guide_legend(title = "Control"))

ggplot(nk_pie, aes(x = "", y = nk_pie, fill = group, alpha=density)) + scale_alpha_continuous(range=c(0,1),limits=c(0,1))+
  geom_col() +
  coord_polar(theta = "y")+ theme_void()+  guides(fill = guide_legend(title = "Treatment"))
dev.off()


# Figure 3C
pdf("./cd8_violinplot_final.pdf",width=30,height=30)
features <- c("Pdcd1","Lag3")
VlnPlot(cd8.combined,features=features,ncol=1,pt.size=0)
features <- c("Ctla4","Havcr2")
VlnPlot(cd8.combined,features=features,ncol=1,pt.size=0)
features <- c("Eomes","Gzmb")
VlnPlot(cd8.combined,features=features,ncol=1,pt.size=0)
features <- c("Gzmk","Tcf7")
VlnPlot(cd8.combined,features=features,ncol=1,pt.size=0)
features <- c("Sell","Bcl2")
VlnPlot(cd8.combined,features=features,ncol=1,pt.size=0)
features <- c("Cd69","Tox")
VlnPlot(cd8.combined,features=features,ncol=1,pt.size=0)
features <- c("Cd44","Klrb1b")
VlnPlot(cd8.combined,features=features,ncol=1,pt.size=0)
dev.off()


# Comment 1
pdf("./featureplot_cd8_comment1_new.pdf")
DefaultAssay(cd8.combined) <- "RNA"
DefaultAssay(cd8.combined.iso) <- "RNA"
DefaultAssay(cd8.combined.nk) <- "RNA"
FeaturePlot(cd8.combined, features=c("Klrk1", "Klrb1b"),split.by="orig.ident")
FeaturePlot(cd8.combined, features=c("Klrb1c", "Hcst"),split.by="orig.ident")
FeaturePlot(cd8.combined, features=c("Tyrobp", "Ncr1"),split.by="orig.ident")
FeaturePlot(cd8.combined, features=c("Tnfsf10", "Tnfrsf1a"),split.by="orig.ident")

FeaturePlot(cd8.combined.iso, features=c("Klrk1", "Klrb1b"))
FeaturePlot(cd8.combined.iso, features=c("Klrb1c", "Hcst"))
FeaturePlot(cd8.combined.iso, features=c("Tyrobp", "Ncr1"))
FeaturePlot(cd8.combined.iso, features=c("Tnfsf10", "Tnfrsf1a"))

FeaturePlot(cd8.combined.nk, features=c("Klrk1", "Klrb1b"))
FeaturePlot(cd8.combined.nk, features=c("Klrb1c", "Hcst"))
FeaturePlot(cd8.combined.nk, features=c("Tyrobp", "Ncr1"))
FeaturePlot(cd8.combined.nk, features=c("Tnfsf10", "Tnfrsf1a"))

dev.off()

# Comment 5
pdf("./featureplot_cd8_comment5.pdf")
DefaultAssay(cd8.combined) <- "RNA"
DefaultAssay(cd8.combined.iso) <- "RNA"
DefaultAssay(cd8.combined.nk) <- "RNA"
FeaturePlot(cd8.combined, features=c("Klrg1", "Cxcr6"))
FeaturePlot(cd8.combined, features=c("Cxcr3", "Itgae"))
FeaturePlot(cd8.combined, features=c("Itga1", "S1pr1"))
FeaturePlot(cd8.combined, features=c("Junb", "Ccr7"))

FeaturePlot(cd8.combined.iso, features=c("Klrg1", "Cxcr6"))
FeaturePlot(cd8.combined.iso, features=c("Cxcr3", "Itgae"))
FeaturePlot(cd8.combined.iso, features=c("Itga1", "S1pr1"))
FeaturePlot(cd8.combined.iso, features=c("Junb", "Ccr7"))

FeaturePlot(cd8.combined.nk, features=c("Klrg1", "Cxcr6"))
FeaturePlot(cd8.combined.nk, features=c("Cxcr3", "Itgae"))
FeaturePlot(cd8.combined.nk, features=c("Itga1", "S1pr1"))
FeaturePlot(cd8.combined.nk, features=c("Junb", "Ccr7"))
dev.off()



#Supp Figure S8 - heatmap
gene_list_cl5nk <- c('Bcl2','Tmsb10','Pdcd1','Ly6c2','Havcr2','H2afz','Ptms','Lag3','Ms4a4c','Gas5',
                     'Ctla2a','Ifng','Icos','Rpl36','Hmgb2','Rpl18a','Cox17','Ccl5','Jund','Rps11',
                     '2410006H16Rik','Pim1','Rgs16','Ly6e','Lmnb1','Cenpa','Tuba1b','Tnfsf8','Bcl2a1d','Rpl37')

gene_list_cl0 <- c('Gzmf','Gzmc','AA467197','Rpl41','Ifitm2','Tpt1','Bcl2a1d','Lgals1','Ifitm1','Aldoa','Ifitm3','Hist1h1b','H2-T23','mt-Co2','Ndufa4','Pfn1','Tnfaip3','mt-Atp6','Pkm','Nt5e','Tox','Vim','Gm42418','mt-Nd3','mt-Co3','Gzmd','Bcl2a1b','Erh','Ctsd','Jund')
gene_list_cl2 <-c('Rpl41',"AA467197","Ifitm2","Gzmf","Ccl5","mt-Nd3",'H2-T23','Gzmc','Rpl36','Tox','Itga4','mt-Co3','mt-Atp6','Bcl2a1d','Ndufa4','mt-Co2','Hmgn2','Rps8','Pfn1','Cd8b1','Ifitm1','Gm8369','Erh',"Gzmk",'Bcl2a1b','Ctla2a','Rbm3',"Ubb","Lgals1","mt-Nd2")
gene_list_isonk<- c('Itga4','Tox','Ifitm2','Rpl41','Ccl5','Gzmf','Gzmc','mt-Atp6','mt-Co3','mt-Nd3','AA467197','H2-T23','Gzmk','Ifitm1','Ifi203','mt-Co2','Gm8369','Ifit3','Ms4a4b','Samhd1','Lgals1','9930111J21Rik2','Prkca','Hist1h2ap','Ifit3b','mt-Cytb','mt-Nd2','Ube2l6','Ifi213','mt-Nd1')

Idents(cd8.combined) <- "seurat_clusters"


pdf("./heatmap_new.pdf")
DoHeatmap(cd8.combined.nk,features = gene_list_cl5nk,group.by="seurat_clusters") + ggtitle("5nk vs rest of nk") 

DoHeatmap(cd8.combined_0,features = gene_list_cl0,group.by="clus.ident") + ggtitle("0nk vs 0iso")
DoHeatmap(cd8.combined_2,features = gene_list_cl2,group.by="clus.ident")+ ggtitle("2nk vs 2iso")
DoHeatmap(cd8.combined,features = gene_list_isonk,group.by="orig.ident")+ ggtitle("nk vs iso")

dev.off()


#Supp Figure S9 - Volcano plot
Idents(cd8.combined) <- "seurat_clusters"
cd8.combined_2 <- subset(cd8.combined, idents = "2",invert = FALSE)
cd8.combined_0 <- subset(cd8.combined, idents = "0",invert = FALSE)

DefaultAssay(cd8.combined.nk) <- "RNA"
unique(cd8.combined.nk$seurat_clusters)

#5nk vs rest
Idents(cd8.combined) <- "clus.ident"
Idents(cd8.combined.nk) <- "clus.ident"
nk.marker <- FindMarkers(cd8.combined.nk, ident.1 = "5_cd8nk", logfc.threshold = 0, 
                         min.cells.group = 1,
                         min.cells.feature = 1,
                         min.pct = 0,verbose = FALSE)
write.csv(nk.marker,paste0("./marker5","_15_cd8nk_all.csv"), row.names = TRUE)
nk.marker <-fread("./marker5_15_cd8nk_all.csv",header=TRUE)
gene_list_cl5nk <- nk.marker$V1[1:30]


nk.marker_cl5_all <- nk.marker %>%
  mutate(pol = ifelse(avg_log2FC > 0, "Down-regulated","Up-regulated"),
         sig = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25, TRUE, FALSE))
nk.marker_cl5_all$X <- rownames(nk.marker_cl5_all)

label_cl5_all <-c("Bcl2","Pdcd1","Ly6c2","Havcr2","Lag3","Ctla2a","Ifng","Icos","Ccl5","Jund","Pim1","Ly6e","Tnfsf8","Bcl2a1d","Satb1","Prf1","Pdcd4","Ctla4","Tox")
v_cl5_all <- ggplot(data = nk.marker_cl5_all,
                    aes(x = -(avg_log2FC),
                        y = -log10(p_val_adj + .Machine$double.xmin),
                        color = sig)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) + 
  geom_label_repel(aes(label = ifelse(X %in% label_cl5_all,as.character(X),""),
                       color = pol),
                   box.padding   = 0.50, 
                   point.padding = 0.30,
                   label.size = 0.20) + 
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_y_continuous(expand = c(0.01,0.01)) +
  xlab("Average log-fold change") + 
  ylab("-log10(Adjusted p-value)") + 
  theme(axis.text = element_text(color = "black",size = 12))
v_cl5_all

#iso vs nk
Idents(cd8.combined) <- "orig.ident"
nk.iso.marker_all <- FindMarkers(cd8.combined, ident.1 = "cd8iso", ident.2 = "cd8nk",logfc.threshold=0,
                                 min.cells.group = 1,
                                 min.cells.feature = 1,pseudocount.use = 0.01,
                                 min.pct = 0,verbose = FALSE)
write.csv(nk.iso.marker_all,"cd8_iso_nk.csv", row.names = TRUE)
nk.iso.marker_all <-fread("./cd8_iso_nk.csv",header=TRUE)
gene_list_isonk <- nk.iso.marker_all$V1[1:30]

nkiso <- nk.iso.marker_all %>%
  mutate(pol = ifelse(avg_log2FC > 0, "Down-regulated","Up-regulated"),
         sig = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25, TRUE, FALSE))
nkiso$X <- rownames(nkiso)
top_up <- nkiso %>%
  filter(avg_log2FC < 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
top_down <- nkiso %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)

label_all <- c(top_down,top_up)
label_all <-c("Itga4","Tox","Ifitm2","Rpl41","Ccl5","Gzmf","Gzmc","mt-Atp6","mt-Co3","mt-Nd3","Gzmk","Ifitm1","Ifi203","Prkca","Bst2","Btg1","Cxcr3","Gzmd","Cd28","Jak2")
v_all <- ggplot(data = nkiso,
                aes(x = -(avg_log2FC),
                    y = -log10(p_val_adj + .Machine$double.xmin),
                    color = sig)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) + 
  geom_label_repel(aes(label = ifelse(X %in% label_all,as.character(X),""),
                       color = pol),
                   box.padding   = 0.50, 
                   point.padding = 0.30,
                   label.size = 0.20, max.overlaps = 1000) +
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_y_continuous(expand = c(0.01,0.01)) +
  xlab("Average log-fold change") + 
  ylab("-log10(Adjusted p-value)") + 
  theme(axis.text = element_text(color = "black",size = 12))
v_all


#0 iso vs 0 nk
nk.iso.marker <- FindMarkers(cd8.combined, ident.1="0_cd8iso", ident.2="0_cd8nk",logfc.threshold=0,
                             min.cells.group = 1,
                             min.cells.feature = 1,pseudocount.use = 0.01,
                             min.pct = 0,verbose = FALSE)

write.csv(nk.iso.marker,"cd8cl0_nk_iso.csv", row.names = TRUE)
nk.iso.marker <-fread("./cd8cl0_nk_iso.csv",header=TRUE)
gene_list_cl0 <- nk.iso.marker$V1[1:30]

nk.marker_cl0 <- nk.iso.marker %>%
  mutate(pol = ifelse(avg_log2FC > 0, "Down-regulated","Up-regulated"),
         sig = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25, TRUE, FALSE))
nk.marker_cl0$X <- rownames(nk.marker_cl0)
label_cl0_all <-c("Gzmf","Gzmc","Ifitm2","Bcl2a1d","Ifitm1","Ifitm3","Pkm","Tox","Vim","Gzmd","Jund","Gzmb","Junb","Itga4","Ly6c2","Fosb","Gzme","Gzma","Prf1","Il2ra","Bcl2")
v_cl0_all <- ggplot(data = nk.marker_cl0,
                    aes(x = -(avg_log2FC),
                        y = -log10(p_val_adj + .Machine$double.xmin),
                        color = sig)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) + 
  geom_label_repel(aes(label = ifelse(X %in% label_cl0_all,as.character(X),""),
                       color = pol),
                   box.padding   = 0.50, 
                   point.padding = 0.30,
                   label.size = 0.20, max.overlaps = 1000) +
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_y_continuous(expand = c(0.01,0.01)) +
  xlab("Average log-fold change") + 
  ylab("-log10(Adjusted p-value)") + 
  theme(axis.text = element_text(color = "black",size = 12))
v_cl0_all

#2 iso vs 2 nk
cluster.ident <- FindMarkers(cd8.combined, ident.1 = "2_cd8iso", ident.2 = "2_cd8nk",logfc.threshold=0,
                             min.cells.group = 1,
                             min.cells.feature = 1,pseudocount.use = 0.01,
                             min.pct = 0,verbose = FALSE)
write.csv(cluster.ident,"cd8cl2_nk_iso.csv", row.names = TRUE)
cluster.ident <-fread("./cd8cl2_nk_iso.csv",header=TRUE)
gene_list_cl2 <- cluster.ident$V1[1:30]

nkiso <- cluster.ident %>%
  mutate(pol = ifelse(avg_log2FC > 0, "Down-regulated","Up-regulated"),
         sig = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25, TRUE, FALSE))
nkiso$X <- rownames(nkiso)
top_up <- nkiso %>%
  filter(avg_log2FC < 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)
top_down <- nkiso %>%
  filter(avg_log2FC > 0) %>%
  top_n(n = 5, wt = -log10(p_val_adj + .Machine$double.xmin)) %>%
  pull(X)

label_all <- c(top_down,top_up)
v2_all <- ggplot(data = nkiso,
                 aes(x = -(avg_log2FC),
                     y = -log10(p_val_adj + .Machine$double.xmin),
                     color = sig)) + 
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("#F8766D","grey65","black","#00BFC4")) + 
  geom_label_repel(aes(label = ifelse(X %in% label_all,as.character(X),""),
                       color = pol),
                   box.padding   = 0.50, 
                   point.padding = 0.30,
                   label.size = 0.20, max.overlaps = 1000) +
  theme_classic() + 
  theme(legend.position = "none") + 
  scale_y_continuous(expand = c(0.01,0.01)) +
  xlab("Average log-fold change") + 
  ylab("-log10(Adjusted p-value)") + 
  theme(axis.text = element_text(color = "black",size = 12))
v2_all

options(future.globals.maxSize= 891289600)


#slingshot
load("afternorm15_cd8inte_upperbound.RData")
Idents(cd8.combined) <- "seurat_clusters"

#removing 10 and 11 cluster
cd8.combined <- subset(cd8.combined, idents = c(10,11), invert = TRUE)
cd8.combined <- subset(cd8.combined, idents = c(7,9), invert = TRUE)

# load into slingshot
p1 <- DimPlot(cd8.combined, reduction = "umap",label=TRUE)
g <- ggplot_build(p1)
color_list <- as.vector(unique(g$data[[1]]["colour"]))

umap <- cd8.combined@reductions$umap@cell.embeddings
seurat_clusters <- cd8.combined@meta.data$seurat_clusters

# run slingshot

fit <- slingshot( data=umap, clusterLabels=seurat_clusters )

pseudotime <- slingPseudotime(fit)


# visualization: clusters
plot( umap, col=seurat_clusters, pch=18, cex=0.7 )
lines( SlingshotDataSet(fit), type = 'lineages' )



seruat_clusters<-as.data.frame(seurat_clusters)
seruat_clusters$color <- color_list[2,]
seruat_clusters$color[which(seruat_clusters$seurat_clusters == 1)] <- color_list[1,]
seruat_clusters$color[which(seruat_clusters$seurat_clusters == 2)] <- color_list[4,]
seruat_clusters$color[which(seruat_clusters$seurat_clusters == 3)] <- color_list[9,]
seruat_clusters$color[which(seruat_clusters$seurat_clusters == 4)] <- color_list[8,]
seruat_clusters$color[which(seruat_clusters$seurat_clusters == 5)] <- color_list[3,]
seruat_clusters$color[which(seruat_clusters$seurat_clusters == 6)] <- color_list[5,]
seruat_clusters$color[which(seruat_clusters$seurat_clusters == 8)] <- color_list[7,]

# Figure 3D
pdf("./slingshot_final.pdf")
plot( umap, col=seruat_clusters$color, pch=18, cex=0.7 )
lines( SlingshotDataSet(fit), type = 'lineages' )


plot( umap, col=seruat_clusters$color, pch=18, cex=0.7 )
lines( SlingshotDataSet(fit), lwd=5) 
dev.off()
# visualization: pseudotime

colors <- colorRampPalette(c("lightyellow", "orange", "darkred"))(100)
plotcol <- colors[cut(pseudotime, breaks=100)]

plot( umap, col=plotcol, pch=18, cex=0.7 )
lines( SlingshotDataSet(fit), lwd=10 )

# visualization: marker expression

DefaultAssay(cd8.combined) <- "RNA"
scale.data <- as.matrix(GetAssayData(cd8.combined, slot = "data"))

list.genes <- c( 
  "Gzmk", "Tox","Pdcd1","Ly6c2")
which(rownames(scale.data) == list.genes[2])
pdf("./slingshot_featureplot.pdf")

for ( i in 1:length(list.genes) ) {
  
  data.i <- scale.data[ rownames(scale.data) == list.genes[i], 
                        match( rownames(umap), colnames(scale.data) ) ]
  
  colors <- colorRampPalette(c("lightyellow", "orange", "darkred"))(100)
  plotcol <- colors[cut(data.i, breaks=100)]
  
  plot( umap, col=plotcol, pch=18, cex=0.7, main=list.genes[i] )
  lines( SlingshotDataSet(fit), lwd=10 )
}
dev.off()


#############################################################################################
#NRAS

library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(igraph)
library(readr)
library(data.table)
library(cowplot)
library(readxl)

nras <- readRDS("~/cellranger/run_cellranger_count_NKcell/NRAS/subset_integrated_Immune.rds")
table(nras$orig.ident)
table(nras$Sample_type)
DimPlot(nras, reduction = "umap.scvi", group.by = "seurat_clusters",label=TRUE)
DefaultAssay(nras) <- "RNA"
nras <- NormalizeData(nras, normalization.method = "LogNormalize", scale.factor = 10000)
nras <- FindVariableFeatures(nras, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
nras <- ScaleData(nras, features = rownames(nras))
saveRDS(nras,file="~/cellranger/run_cellranger_count_NKcell/NRAS/subset_integrated_Immune_normalize.rds")
DefaultAssay(nras) <- "SCT"
nras <- ScaleData(nras, features = rownames(nras))
#sctype usage
# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = nras[["RNA"]]$scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
save(es.max, file="~/cellranger/run_cellranger_count_NKcell/NRAS/sctype_score.RData")
# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(nras@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(nras@meta.data[nras@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(nras@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
sctype_anno <- sctype_scores[order(sctype_scores$cluster),][1:2]

nras$sctype_anno <- 1
cln<-as.numeric(as.vector(unlist(sctype_anno[,1])))
#i<-1
for(i in cln){
  nras$sctype_anno[which(nras$seurat_clusters==i)] <- as.character(sctype_anno[i+1,2])  
}
table(nras$sctype_anno)
Idents(nras)<-"sctype_anno"
p1 <- DimPlot(nras, reduction = "umap.scvi", label = TRUE, repel = TRUE) + ggtitle('sctype anno')
print(p1)
save(nras,file="~/cellranger/run_cellranger_count_NKcell/NRAS/nras_analysis_sctype.RData")
load("~/cellranger/run_cellranger_count_NKcell/NRAS/nras_analysis_sctype.RData")


Idents(nras) <- "Sample_type"
cd8.ext <- subset(nras, idents = unique(nras$Sample_type)[5:8])
DimPlot(cd8.ext, reduction = "umap.scvi", label = TRUE, repel = TRUE) + ggtitle('Sample type')

Idents(cd8.ext) <- "sctype_anno"
cd8.ext_nras <- subset(cd8.ext, idents = "Memory CD8+ T cells")
DimPlot(cd8.ext_nras, reduction = "umap.scvi", label = TRUE, repel = TRUE) + ggtitle('CD8 only')

save(cd8.ext_nras,file="cd8_nras_analysis_sctype.RData")
load("cd8_nras_analysis_sctype.RData")

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(cd8.ext_nras) <- "RNA"
cd8.ext_nras <- RunPCA(cd8.ext_nras, npcs = 30, verbose = FALSE)
cd8.ext_nras <- RunUMAP(cd8.ext_nras, reduction = "pca", dims = 1:30)
cd8.ext_nras <- FindNeighbors(cd8.ext_nras, reduction = "pca", dims = 1:30)
cd8.ext_nras <- FindClusters(cd8.ext_nras, resolution = 1.1)


#9 clusters
save(cd8.ext_nras,file="after_cd8_cluster_res1.1nras.RData")
load("after_cd8_cluster_res1.1_nras.RData")

# Visualization
# Figure 3F
pdf("./cd8_9clu_nras.pdf")

Idents(cd8.ext_nras) <- "seurat_clusters"
p1 <- DimPlot(cd8.ext_nras, reduction = "umap",label=TRUE,split.by = "Sample_type")
p2 <- DimPlot(cd8.ext_nras, reduction = "umap", label = TRUE, repel = TRUE)
p1
p2


# Figure 3G
preprop_table <- table(cd8.ext_nras$seurat_clusters,cd8.ext_nras$Sample_type)
preprop_table
preprop_table[,1] <- preprop_table[,1]/sum(preprop_table[,1])
preprop_table[,2] <- preprop_table[,2]/sum(preprop_table[,2])
preprop_table[,3] <- preprop_table[,3]/sum(preprop_table[,3])
preprop_table[,4] <- preprop_table[,4]/sum(preprop_table[,4])
preprop_table
prop <- as.data.frame(preprop_table)
library(tidyr)

p3 <- ggplot(prop, aes(x = Var2, y = Freq, fill = as.factor(Var1))) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Each cluster in Each sample",
       x = "Sample type",
       y = "Proportion") +
  theme_minimal() +
  scale_fill_discrete(name = "Row")
print(p2)

print(p1)
print(p3)
dev.off()

# Figure 3F
pdf("./cd8_feature_plot_nras.pdf")
FeaturePlot(cd8.ext_nras, features=c("Havcr2","Bcl2", "Eomes","Gzmb"))
dev.off()

table(cd8.ext_nras$Sample_type)


#Fig 3H

library(ggpubr)
library("readxl")
library(ggplot2)
### cd8.ext_nras is the seurat object
load("after_cd8_cluster_res1.1.RData")

count <- cd8.ext_nras@assays$RNA@layers$counts
UMAP <- cd8.ext_nras@reductions$umap@cell.embeddings


##### NKRT signature_option2: use the gene list provided by Dr.Song
my_data <- read_excel("./NRAS/signature_gene_list.xlsx")

my_data.count <- count[which(toupper(rownames(cd8.ext_nras))%in% toupper(my_data$gene)),]
my_data.count <- t(my_data.count)

my_data <- my_data[which(toupper(my_data$gene)%in% toupper(rownames(cd8.ext_nras))),]
direction <- my_data[,2]
direction <- direction %>%
  mutate(exp = ifelse(exp == "High", 1, -1))%>% 
  as.matrix()

my_data.score <- my_data.count%*%direction

my_data.score.mod <- my_data.score
my_data.score.mod[which(my_data.score.mod> quantile(my_data.score,0.90))] <- quantile(my_data.score,0.90)
my_data.score.mod[which(my_data.score.mod< quantile(my_data.score,0.10))] <- quantile(my_data.score,0.10)


Idents(cd8.ext_nras) <- "Sample_type"
unique(cd8.ext_nras$Sample_type)
con_seu <- subset(cd8.ext_nras, idents=c("NRAScontrol"),invert=FALSE)
pd1_seu <- subset(cd8.ext_nras, idents=c("NRASaPD1"),invert=FALSE)
nk11_seu <- subset(cd8.ext_nras, idents=c("NRASaNK11"),invert=FALSE)
pd1nk11_seu <- subset(cd8.ext_nras, idents=c("NRASaPD1aNK11"),invert=FALSE)

con_count <- con_seu@assays$RNA@layers$counts
con_UMAP <- con_seu@reductions$umap@cell.embeddings

pd1_count <- pd1_seu@assays$RNA@layers$counts
pd1_UMAP <- pd1_seu@reductions$umap@cell.embeddings

nk11_count <- nk11_seu@assays$RNA@layers$counts
nk11_UMAP <- nk11_seu@reductions$umap@cell.embeddings

pd1nk11_count <- pd1nk11_seu@assays$RNA@layers$counts
pd1nk11_UMAP <- pd1nk11_seu@reductions$umap@cell.embeddings

##### NKRT signature_option2: use the gene list provided by Dr.Song
my_data <- read_excel("./NRAS/signature_gene_list.xlsx")

my_data.count_con <- con_count[which(toupper(rownames(con_seu))%in% toupper(my_data$gene)),]
my_data.count_con <- t(my_data.count_con)

my_data.count_pd1 <- pd1_count[which(toupper(rownames(pd1_seu))%in% toupper(my_data$gene)),]
my_data.count_pd1 <- t(my_data.count_pd1)

my_data.count_nk11 <- nk11_count[which(toupper(rownames(nk11_seu))%in% toupper(my_data$gene)),]
my_data.count_nk11 <- t(my_data.count_nk11)

my_data.count_pd1nk11 <- pd1nk11_count[which(toupper(rownames(pd1nk11_seu))%in% toupper(my_data$gene)),]
my_data.count_pd1nk11 <- t(my_data.count_pd1nk11)

my_data <- my_data[which(toupper(my_data$gene)%in% toupper(rownames(pd1nk11_seu))),]
direction <- my_data[,2]
direction <- direction %>%
  mutate(exp = ifelse(exp == "High", 1, -1))%>% 
  as.matrix()

my_data.score_con <- my_data.count_con%*%direction
my_data.score_pd1 <- my_data.count_pd1%*%direction
my_data.score_nk11 <- my_data.count_nk11%*%direction
my_data.score_pd1nk11 <- my_data.count_pd1nk11%*%direction


my_data.score.mod_con <- my_data.score_con
my_data.score.mod_con[which(my_data.score.mod_con> quantile(my_data.score_con,0.90))] <- quantile(my_data.score_con,0.90)
my_data.score.mod_con[which(my_data.score.mod_con< quantile(my_data.score_con,0.10))] <- quantile(my_data.score_con,0.10)


my_data.score.mod_pd1 <- my_data.score_pd1
my_data.score.mod_pd1[which(my_data.score.mod_pd1> quantile(my_data.score_pd1,0.90))] <- quantile(my_data.score_pd1,0.90)
my_data.score.mod_pd1[which(my_data.score.mod_pd1< quantile(my_data.score_pd1,0.10))] <- quantile(my_data.score_pd1,0.10)


my_data.score.mod_nk11 <- my_data.score_nk11
my_data.score.mod_nk11[which(my_data.score.mod_nk11> quantile(my_data.score_nk11,0.90))] <- quantile(my_data.score_nk11,0.90)
my_data.score.mod_nk11[which(my_data.score.mod_nk11< quantile(my_data.score_nk11,0.10))] <- quantile(my_data.score_nk11,0.10)


my_data.score.mod_pd1nk11 <- my_data.score_pd1nk11
my_data.score.mod_pd1nk11[which(my_data.score.mod_pd1nk11> quantile(my_data.score_pd1nk11,0.90))] <- quantile(my_data.score_pd1nk11,0.90)
my_data.score.mod_pd1nk11[which(my_data.score.mod_pd1nk11< quantile(my_data.score_pd1nk11,0.10))] <- quantile(my_data.score_pd1nk11,0.10)


####Violin plot
pdf("./NKRT_signature_violin_with_dots.pdf")
changed <- my_data.score_con[,1]+(-1*summary(my_data.score_con[,1])[1])
df <- data.frame(changed)
df$condition <- "Control"


changed <- my_data.score_pd1[,1]+(-1*summary(my_data.score_pd1[,1])[1])
df1 <- data.frame(changed)
df1$condition <- "PD1"


changed <- my_data.score_nk11[,1]+(-1*summary(my_data.score_nk11[,1])[1])
df2 <- data.frame(changed)
df2$condition <- "NK11"


changed <- my_data.score_pd1nk11[,1]+(-1*summary(my_data.score_pd1nk11[,1])[1])
df3 <- data.frame(changed)
df3$condition <- "PD1NK11"
df <- rbind(df,df1,df2,df3)
ggplot(df, aes(x=condition, y=changed,fill=condition)) + 
  geom_violin(position = position_dodge(width = 3)) +
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 1))+ ggtitle("Raw")+ stat_compare_means(method = "kruskal.test", label.y = 9)


changed <- log(my_data.score_con[,1]+(-1*summary(my_data.score_con[,1])[1])+1)
df <- data.frame(changed)
df$condition <- "Control"



changed <- log(my_data.score_pd1[,1]+(-1*summary(my_data.score_pd1[,1])[1])+1)
df1 <- data.frame(changed)
df1$condition <- "PD1"

changed <- log(my_data.score_nk11[,1]+(-1*summary(my_data.score_nk11[,1])[1])+1)
df2 <- data.frame(changed)
df2$condition <- "NK11"


changed <- log(my_data.score_pd1nk11[,1]+(-1*summary(my_data.score_pd1nk11[,1])[1])+1)
df3 <- data.frame(changed)
df3$condition <- "PD1NK11"
df <- rbind(df,df1,df2,df3)

ggplot(df, aes(x=condition, y=changed,fill=condition)) + 
  geom_violin(position = position_dodge(width = 3)) +
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 1))+ ggtitle("log transform")


changed <- my_data.score.mod_con[,1]+(-1*summary(my_data.score.mod_con[,1])[1])
df <- data.frame(changed)
df$condition <- "Control"


changed <- my_data.score.mod_pd1[,1]+(-1*summary(my_data.score.mod_pd1[,1])[1])
df1 <- data.frame(changed)
df1$condition <- "PD1"

changed <- my_data.score.mod_nk11[,1]+(-1*summary(my_data.score.mod_nk11[,1])[1])
df2 <- data.frame(changed)
df2$condition <- "NK11"


changed <- my_data.score.mod_pd1nk11[,1]+(-1*summary(my_data.score.mod_pd1nk11[,1])[1])
df3 <- data.frame(changed)
df3$condition <- "PD1NK11"
df <- rbind(df,df1,df2,df3)

ggplot(df, aes(x=condition, y=changed,fill=condition)) + 
  geom_violin(position = position_dodge(width = 3)) +
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 1))+ ggtitle("raw - top, bottom 10% quantile adjusted")

changed <- log(my_data.score.mod_con[,1]+(-1*summary(my_data.score.mod_con[,1])[1])+1)
df <- data.frame(changed)
df$condition <- "Control"


changed <- log(my_data.score.mod_pd1[,1]+(-1*summary(my_data.score.mod_pd1[,1])[1])+1)
df1 <- data.frame(changed)
df1$condition <- "PD1"


changed <- log(my_data.score.mod_nk11[,1]+(-1*summary(my_data.score.mod_nk11[,1])[1])+1)
df2 <- data.frame(changed)
df2$condition <- "NK11"

changed <- log(my_data.score.mod_pd1nk11[,1]+(-1*summary(my_data.score.mod_pd1nk11[,1])[1])+1)
df3 <- data.frame(changed)
df3$condition <- "PD1NK11"
df <- rbind(df,df1,df2,df3)

ggplot(df, aes(x=condition, y=changed,fill=condition)) + 
  geom_violin(position = position_dodge(width = 3)) +
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 1))+ ggtitle("log transform - top, bottom 10% quantile adjusted")
dev.off()



pdf("./NKRT_signature_violin.pdf")
changed <- my_data.score_con[,1]+(-1*summary(my_data.score_con[,1])[1])
df <- data.frame(changed)
df$condition <- "Control"


changed <- my_data.score_pd1[,1]+(-1*summary(my_data.score_pd1[,1])[1])
df1 <- data.frame(changed)
df1$condition <- "PD1"


changed <- my_data.score_nk11[,1]+(-1*summary(my_data.score_nk11[,1])[1])
df2 <- data.frame(changed)
df2$condition <- "NK11"


changed <- my_data.score_pd1nk11[,1]+(-1*summary(my_data.score_pd1nk11[,1])[1])
df3 <- data.frame(changed)
df3$condition <- "PD1NK11"
df <- rbind(df,df1,df2,df3)
ggplot(df, aes(x=condition, y=changed,fill=condition)) + 
  geom_violin()+ ggtitle("Raw")+ stat_compare_means(method = "kruskal.test", label.y = 400)


changed <- log(my_data.score_con[,1]+(-1*summary(my_data.score_con[,1])[1])+1)
df <- data.frame(changed)
df$condition <- "Control"



changed <- log(my_data.score_pd1[,1]+(-1*summary(my_data.score_pd1[,1])[1])+1)
df1 <- data.frame(changed)
df1$condition <- "PD1"

changed <- log(my_data.score_nk11[,1]+(-1*summary(my_data.score_nk11[,1])[1])+1)
df2 <- data.frame(changed)
df2$condition <- "NK11"


changed <- log(my_data.score_pd1nk11[,1]+(-1*summary(my_data.score_pd1nk11[,1])[1])+1)
df3 <- data.frame(changed)
df3$condition <- "PD1NK11"
df <- rbind(df,df1,df2,df3)

ggplot(df, aes(x=condition, y=changed,fill=condition)) + 
  geom_violin()+ ggtitle("log transform")+ stat_compare_means(method = "kruskal.test", label.y = 400)



changed <- my_data.score.mod_con[,1]+(-1*summary(my_data.score.mod_con[,1])[1])
df <- data.frame(changed)
df$condition <- "Control"


changed <- my_data.score.mod_pd1[,1]+(-1*summary(my_data.score.mod_pd1[,1])[1])
df1 <- data.frame(changed)
df1$condition <- "PD1"

changed <- my_data.score.mod_nk11[,1]+(-1*summary(my_data.score.mod_nk11[,1])[1])
df2 <- data.frame(changed)
df2$condition <- "NK11"


changed <- my_data.score.mod_pd1nk11[,1]+(-1*summary(my_data.score.mod_pd1nk11[,1])[1])
df3 <- data.frame(changed)
df3$condition <- "PD1NK11"
df <- rbind(df,df1,df2,df3)

ggplot(df, aes(x=condition, y=changed,fill=condition)) + 
  geom_violin()+ ggtitle("raw - top, bottom 10% quantile adjusted")+ stat_compare_means(method = "kruskal.test", label.y = 400)


changed <- log(my_data.score.mod_con[,1]+(-1*summary(my_data.score.mod_con[,1])[1])+1)
df <- data.frame(changed)
df$condition <- "Control"


changed <- log(my_data.score.mod_pd1[,1]+(-1*summary(my_data.score.mod_pd1[,1])[1])+1)
df1 <- data.frame(changed)
df1$condition <- "PD1"


changed <- log(my_data.score.mod_nk11[,1]+(-1*summary(my_data.score.mod_nk11[,1])[1])+1)
df2 <- data.frame(changed)
df2$condition <- "NK11"

changed <- log(my_data.score.mod_pd1nk11[,1]+(-1*summary(my_data.score.mod_pd1nk11[,1])[1])+1)
df3 <- data.frame(changed)
df3$condition <- "PD1NK11"
df <- rbind(df,df1,df2,df3)

ggplot(df, aes(x=condition, y=changed,fill=condition)) + 
  geom_violin()+ ggtitle("log transform - top, bottom 10% quantile adjusted")+ stat_compare_means(method = "kruskal.test", label.y = 400)

dev.off()



pdf("./NKRT_signature_violin_wilcoxon.pdf")
changed <- my_data.score_con[,1]+(-1*summary(my_data.score_con[,1])[1])
df <- data.frame(changed)
df$condition <- "Control"


changed <- my_data.score_pd1[,1]+(-1*summary(my_data.score_pd1[,1])[1])
df1 <- data.frame(changed)
df1$condition <- "PD1"


changed <- my_data.score_nk11[,1]+(-1*summary(my_data.score_nk11[,1])[1])
df2 <- data.frame(changed)
df2$condition <- "NK11"


changed <- my_data.score_pd1nk11[,1]+(-1*summary(my_data.score_pd1nk11[,1])[1])
df3 <- data.frame(changed)
df3$condition <- "PD1NK11"
df <- rbind(df,df1,df2,df3)
ggplot(df, aes(x=condition, y=changed,fill=condition)) + 
  geom_violin()+ ggtitle("Raw")+ stat_compare_means(comparisons = list(c("Control", "PD1"), 
                                                                       c("Control", "NK11"), 
                                                                       c("Control", "PD1NK11"), 
                                                                       c("PD1", "NK11"), 
                                                                       c("PD1", "PD1NK11"), 
                                                                       c("NK11", "PD1NK11")),
                                                    method = "wilcox.test", label = "p.signif", 
                                                    p.adjust.method = "bonferroni", 
                                                    label.y = c(250, 300, 450, 270, 470, 500))  # Adjust label.y for better visualization

ggplot(df, aes(x=condition, y=changed,fill=condition)) + 
  geom_violin()+ ggtitle("Raw")+ stat_compare_means(comparisons = list(c("Control", "PD1"), 
                                                                       c("Control", "NK11"), 
                                                                       c("Control", "PD1NK11"), 
                                                                       c("PD1", "NK11"), 
                                                                       c("PD1", "PD1NK11"), 
                                                                       c("NK11", "PD1NK11")),
                                                    method = "wilcox.test", label = "p.format",
                                                    p.adjust.method = "bonferroni", 
                                                    label.y = c(250, 310, 450, 290, 470, 500))  # Adjust label.y for better visualization


dev.off()
#################################################################
######################  Figure 4I: reanalysis of public data (Daniel et al) #############
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Matrix)

DAT <- Read10X(data.dir = '/fs/ess/PAS2579/NoJoon/GSE188666/out',gene.column=1)
DAT <- CreateSeuratObject(counts = DAT,project = "GSE18866",min.cells = 0, min.features = 0)

a <-sapply(strsplit(colnames(DAT),'_',1),'[',3)
DAT$time <- a


meta <- read.table('/fs/ess/PAS2579/NoJoon/GSE188666/out/GSE188666_scRNA_LCMV_metadata.tsv')

identical(colnames(DAT),rownames(meta))  # check if the order of cells are the same

DAT <- AddMetaData(DAT, meta)

Idents(DAT) <- DAT$ident
UMAP <- data.frame(UMAP_1=DAT$UMAP_1,UMAP_2=DAT$UMAP_2)
UMAP <- as.matrix(UMAP)
DAT[['umap']] <-CreateDimReducObject(embeddings = UMAP,key = 'UMAP_',assay = DefaultAssay(DAT))
DAT$clusters <- DAT$ident


DimPlot(DAT)

count <- DAT@assays$RNA@counts
UMAP <- DAT@reductions$umap@cell.embeddings

### NKRT signature
GENE <- read.csv('/fs/ess/PAS2579/NoJoon/geneList2_song.csv')   # alternative set2
GENE.count <- count[which(toupper(rownames(count))%in% toupper(GENE$gene)),]
GENE.count <- t(GENE.count)

a <- rep(1,ncol(GENE.count))
a <- as.matrix(a)
GENE.score <- GENE.count%*%a

summary(GENE.score[,1])  
hist(GENE.score[,1])
hist(log(GENE.score[,1]+1))


rownames(GENE.score) <- rownames(GENE.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(GENE.score)),]

df <- data.frame(UMAP.sub,GENE.score, log(GENE.score+1))

mid <- median(log(GENE.score))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(GENE.score)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()



## RAR

RAR.target <- read.csv('/fs/ess/PAS2579/NoJoon/set3.csv',header=T)

RAR.count <- count[which(rownames(count)%in% RAR.target$Set3_RAR),]
RAR.count <- t(RAR.count)

a <- rep(1,ncol(RAR.count))
a <- as.matrix(a)
RAR.score <- RAR.count%*%a

summary(RAR.score[,1])


rownames(RAR.score) <- rownames(RAR.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(RAR.score)),]

df <- data.frame(UMAP.sub,RAR.score)
mid <- median(log(RAR.score))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(RAR.score)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()



##### combine RAR and NKRT together

RAR <- RAR.target$Set3_RAR


TARGET.GENE <- c(GENE$gene,toupper(RAR))

TARGET.count <- count[which(toupper(rownames(count))%in% TARGET.GENE),]
TARGET.count <- t(TARGET.count)


a <- vector()

IND <- which(colnames(TARGET.count)%in% toupper(GENE$gene))

for(i in 1:ncol(TARGET.count)){
  if (i %in% IND){
    temp <- colnames(TARGET.count)[i]
    tt <- GENE[which(toupper(GENE$gene)==temp),]$exp
    a[i] <- ifelse(tt=='High',0.1,0.1)
  }else{
    a[i] <- 0.1
  }
}



a <- as.matrix(a)
TARGET.score <- TARGET.count %*% a

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(TARGET.score)),]
TARGET.score.mod <- log(TARGET.score[,1])

for (i in 1:length(TARGET.score.mod)){
  TARGET.score.mod[i] <- ifelse(TARGET.score.mod[i]>4.55,4.55,TARGET.score.mod[i])
  TARGET.score.mod[i] <- ifelse(TARGET.score.mod[i]<3.25,3.25,TARGET.score.mod[i])
}

hist(TARGET.score.mod)


df <- data.frame(UMAP.sub,TARGET.score,TARGET.score.mod)

mid <- median(TARGET.score.mod)
ggplot(df,aes(UMAP_1, UMAP_2,color=TARGET.score.mod))+geom_point(size=0.2) +
  scale_colour_gradient2(low = "royalblue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()

#################################################################
################## Figure 6 F,G,H (zheng et al., GSE221064)
## create seurat object
library(Seurat)
library(dplyr)
library(ggplot2)

dir.in <- '/fs/ess/PAS2579/NoJoon/GSE221064/data_seuratV4/'

NAMES <- list.files(dir.in)

data.id <- grep('GSM',NAMES)
DATA.names <- NAMES[data.id]

for (i in 1:length(data.id)){
  data.dir <- paste0(dir.in, DATA.names[i])
  data <- Read10X(data.dir = data.dir)
  seurat.obj <- CreateSeuratObject(counts = data,project = DATA.names[i],min.features = 200,min.cells = 3)
  seurat.obj[['percent.mt']] <- PercentageFeatureSet(seurat.obj,pattern = '^mt-')
  
  saveRDS(seurat.obj,file = paste0(dir.in,'seurat/',DATA.names[i],'_seurat.rds'))
}


dir.in <- '/fs/ess/PAS2579/NoJoon/GSE221064/data_seuratV4/seurat/'
FILES <-list.files(dir.in,pattern = 'rds')
FILES


i <-1
DATA1 <- readRDS(paste0(dir.in,FILES[i]))
DATA1 <- subset(DATA1,subset = percent.mt <20 & nCount_RNA <75000)

i <-2
DATA2 <- readRDS(paste0(dir.in,FILES[i]))
DATA2 <- subset(DATA2,subset = percent.mt <25 & nCount_RNA < 100000 )

i <-3
DATA3 <- readRDS(paste0(dir.in,FILES[i]))
DATA3 <- subset(DATA3,subset = percent.mt <20 & nCount_RNA < 100000 )

#### log normalization and combine

DATA1 <-NormalizeData(DATA1, normalization.method = 'LogNormalize',scale.factor = 10000)
DATA1 <-FindVariableFeatures(DATA1,selection.method='vst',nfeatures=2000)
DATA1 <-RenameCells(DATA1,add.cell.id = 'A7')  # add prefix to the colnames to avoid duplicated cell names during integration

DATA2 <-NormalizeData(DATA2, normalization.method = 'LogNormalize',scale.factor = 10000)
DATA2 <-FindVariableFeatures(DATA2,selection.method='vst',nfeatures=2000)
DATA2 <-RenameCells(DATA2,add.cell.id = 'A7.1')  # add prefix to the colnames to avoid duplicated cell names during integration

DATA3 <-NormalizeData(DATA3, normalization.method = 'LogNormalize',scale.factor = 10000)
DATA3 <-FindVariableFeatures(DATA3,selection.method='vst',nfeatures=2000)
DATA3 <-RenameCells(DATA3,add.cell.id = 'B7')  # add prefix to the colnames to avoid duplicated cell names during integration

DATA.list <- list(DATA1,DATA2,DATA3)

features <- SelectIntegrationFeatures(object.list = DATA.list, nfeatures = 3000)
DATA.anchors <- FindIntegrationAnchors(object.list = DATA.list,
                                       anchor.features = features)
AB.combined <- IntegrateData(anchorset = DATA.anchors)

saveRDS(AB.combined, file = '/fs/ess/PAS2579/NoJoon/GSE221064/data_seuratV4/seurat/combined/AB.combined.rds')

AB.combined <- readRDS(file = '/fs/ess/PAS2579/NoJoon/GSE221064/data_seuratV4/seurat/combined/AB.combined.rds')

AB.combined <- ScaleData(AB.combined)
AB.combined <- RunPCA(AB.combined, npcs = 30)
ElbowPlot(AB.combined)
AB.combined <- RunUMAP(AB.combined, reduction = "pca", dims = 1:30)

AB.combined <- FindNeighbors(AB.combined,dims = 1:30)
AB.combined <- FindClusters(AB.combined,resolution = 0.3)

# AB.combined <- FindClusters(AB.combined,resolution = 0.8)




###annotate the cells using scType

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
#tissue = "Heart" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = AB.combined[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(AB.combined@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(AB.combined@meta.data[AB.combined@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(AB.combined@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


AB.combined@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  AB.combined@meta.data$customclassif[AB.combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


DimPlot(AB.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  


DefaultAssay(AB.combined) <- 'RNA'
Idents(AB.combined) <- AB.combined$customclassif

DATA.NK <- subset(AB.combined,subset = customclassif =='Natural killer  cells')

DATA.NK <- DATA.NK %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData %>%
  RunPCA %>%
  RunUMAP(dims = 1:20)

DATA.NK <- FindNeighbors(DATA.NK,dims = 1:20)

DATA.NK <- FindClusters(DATA.NK,resolution = 0.3) 

DimPlot(DATA.NK, label = T)

DATA.NK.sub <- subset(DATA.NK, subset = seurat_clusters %in% c(0,1,2,3,4))
Idents(DATA.NK.sub) <- DATA.NK.sub$seurat_clusters
DotPlot(DATA.NK.sub,features = c('Cd69','Eomes','Bcl2','Gzmb'))


#### heatmap 
GENES <- c('Il2ra','Il2rb','Il4r','Il10ra','Il17ra','Il12rb1','Il15ra','Ifngr2','Ifngr1','Ifnar1','Ifnar2')
#GENES <- c('Cd69','Itgam','Sell','Gzmb','Bcl2','Eomes','Tbx21')


averaged <- AverageExpression(DATA.NK.sub, assays = 'RNA',group.by = 'seurat_clusters')

PLOT <- averaged[[1]]
idx <- which(rownames(PLOT)%in% GENES)
FIND <- rownames(PLOT)[idx]
setdiff(GENES,FIND)
PLOT.sub <- as.matrix(PLOT[idx,])

temp <- apply(PLOT.sub,1,sum)
PLOT.sub <- PLOT.sub[-which(temp==0),]

pheatmap(PLOT.sub,scale='row')


################## Figure 7

#################################################################
#Figure 7AB
#Feldman
library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(readr)
library(data.table)
library(cowplot)
library(readxl)

# data directory (change on your computer)
# For PBMC Seurat tutorial data, put "filtered_gene_bc_matrices/" within the data_dir
data_dir <- c("../")
# PBMC Seurat tutorial data
a <- readLines("../../GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt")
a[1] <-sub('.','', a[1])
a <- a[-2]
writeLines(a,"./GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt")

feldman.data <- read.table(file = "../GSE120575_patient_ID_single_cells.txt", sep="\t", header=TRUE, row.names = 1)
rownames(feldman.data) <- feldman.data$title
feldman.data <- feldman.data[,4:6]
feldman.data <- feldman.data %>% filter(feldman.data$characteristics..therapy == 'anti-PD1')
for (i in c(1:16291)){
  feldman.data$patient[i] <- str_split(feldman.data$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.,"_",n=2, simplify=TRUE)[i,2]
}
feldman.data$characteristics..patinet.ID..Pre.baseline..Post..on.treatment. <- gsub("Pre_.*", "Pre", feldman.data$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)
feldman.data$characteristics..patinet.ID..Pre.baseline..Post..on.treatment. <- gsub("Post_.*", "Post", feldman.data$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)
table(feldman.data$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)

feldman.data1 <- read.table(file = "../GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt", sep="\t", dec = ".", header=FALSE, skip=1,row.names = 1)
header <- read.table( file = "../GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt", sep="\t", dec = ".", nrows=1 )
colnames(feldman.data1) <- unlist(header)
feldman.data1 <- feldman.data1 %>% dplyr::select(rownames(feldman.data))

# Initialize the Seurat object with the raw (non-normalized data).
feldman <- CreateSeuratObject(counts = feldman.data1, 
                              project = "NKcells_count", 
                              meta.data = feldman.data,
                              min.cells = 3, #low quality genes
                              min.features = 200) #at least 200 features

#save(feldman, file="./feldman.RData")
save(feldman, file="./feldman_all_therapy.RData")

load('../../feldman/feldman_all_therapy.RData')

Idents(feldman) <- "characteristics..therapy"
unique(feldman$characteristics..therapy)

###only pd1
feldman_pd1 <- subset(feldman,idents = c("anti-PD1"), invert = FALSE)
feldman_pd1$treatment_response <- paste(feldman_pd1@meta.data$characteristics..patinet.ID..Pre.baseline..Post..on.treatment., feldman_pd1@meta.data$characteristics..response, sep="_")

# Normalize counts
feldman_pd1 <- NormalizeData(feldman_pd1)

# Identify the most informative genes
feldman_pd1 <- FindVariableFeatures(feldman_pd1, 
                                    selection.method = "vst", 
                                    nfeatures = 2000)

# Scale the genes
feldman_pd1 <- ScaleData(feldman_pd1)
# Compute PCA (50)
feldman_pd1 <- RunPCA(feldman_pd1, features = VariableFeatures(object = feldman_pd1))
#rectangular data

#From here graph data
# Construct cell-cell nearest neighbors graph
feldman_pd1 <- FindNeighbors(feldman_pd1, dims = 1:10)
# Cluster the graph using Louvain
feldman_pd1 <- FindClusters(feldman_pd1, resolution = 0.5)

# # save counts for singleR
counts <- feldman_pd1[["RNA"]]@counts
metadata <- feldman_pd1@meta.data
#save(counts,metadata,file = "./feldmancounts_norm.RData")
#save(feldman, file="./feldmanafternorm.RData")
save(counts,metadata,file = "./feldman_pd1_norm.RData")
save(feldman_pd1, file="../../feldman/feldman_pd1_therapyafternorm.RData")
load("../../feldman/feldman_pd1_therapyafternorm.RData")

# Figure 7AB

#feldman_pd1_cd8 <- subset(feldman_pd1, idents = c(0,2,6), invert = FALSE)
feldman_pd1_cd8 <- RunUMAP(feldman_pd1, dims = 1:10) #only visualization usage. UMAP  nonlinear / using 10 pca features /hyper parameter
DimPlot(feldman_pd1_cd8, reduction= "umap", label=TRUE)

table(feldman_pd1_cd8@meta.data$treatment_response)
#Post_Non-responder     Post_Responder  Pre_Non-responder      Pre_Responder 
#6334                   1524            2604                   1191 

#change the number of cells accordingly
feldman.pretreatment <- subset( feldman_pd1_cd8, 
                                cells=c( 
                                  base::sample( which( feldman_pd1_cd8@meta.data$treatment_response == "Pre_Non-responder" ), 1191, replace=FALSE ),
                                  base::sample( which( feldman_pd1_cd8@meta.data$treatment_response == "Pre_Responder" ), 1191, replace=FALSE ) ) )


DimPlot(feldman.pretreatment,split.by = "treatment_response", reduction= "umap", label=TRUE)

preprop_table <- table(feldman.pretreatment$seurat_clusters, feldman.pretreatment$treatment_response)
preprop_table
preprop_table[,1] <- preprop_table[,1]/sum(preprop_table[,1])
preprop_table[,2] <- preprop_table[,2]/sum(preprop_table[,2])
preprop_table
preprop <- as.data.frame((preprop_table))#.matrix((preprop_table))
which(preprop$Var1=="7")
preprop <- preprop[which(preprop$Var1=="7"),]
preprop

# Figure 7AB

pdf("../feldman/feldman_all_therapyproportion_pd1only.pdf")
ggplot(data=preprop, aes(x=Var1, y=Freq, fill = Var2)) +
  geom_bar(position="dodge",stat="identity", width=0.5) +
  ylab("Proportion") + 
  xlab("Cell type") +
  theme(text = element_text(size=20))+ guides(fill=guide_legend(title="Response"))
dev.off()
prop.test( c(246,190), c(2725,2725))
###only pd1 ends

#all
feldman$treatment_response <- paste(feldman@meta.data$characteristics..patinet.ID..Pre.baseline..Post..on.treatment., feldman@meta.data$characteristics..response, sep="_")

# Normalize counts
feldman <- NormalizeData(feldman)

# Identify the most informative genes
feldman <- FindVariableFeatures(feldman, 
                                selection.method = "vst", 
                                nfeatures = 2000)

# Scale the genes
feldman <- ScaleData(feldman)
# Compute PCA (50)
feldman <- RunPCA(feldman, features = VariableFeatures(object = feldman))
#rectangular data

#From here graph data
# Construct cell-cell nearest neighbors graph
feldman <- FindNeighbors(feldman, dims = 1:10)
# Cluster the graph using Louvain
feldman <- FindClusters(feldman, resolution = 0.5)

# # save counts for singleR
counts <- feldman[["RNA"]]@counts
metadata <- feldman@meta.data
#save(counts,metadata,file = "./feldmancounts_norm.RData")
#save(feldman, file="./feldmanafternorm.RData")
save(counts,metadata,file = "./feldman_all_therapycounts_norm.RData")
save(feldman, file="../../feldman/feldman_all_therapyafternorm.RData")
load("../../feldman/feldman_all_therapyafternorm.RData",verbose=T)

# Figure 7AB

pdf("./feldman_all_therapy_featureplot.pdf")
# Compute UMAP coordinates for visualization
feldman <- RunUMAP(feldman, dims = 1:10) #only visualization usage. UMAP  nonlinear / using 10 pca features /hyper parameter
DimPlot(feldman, reduction= "umap", label=TRUE)

table(feldman@meta.data$treatment_response)
#Post_Non-responder     Post_Responder  Pre_Non-responder      Pre_Responder 
#6334                   1524            2604                   1191 

#All therapy
#Post_Non-responder     Post_Responder  Pre_Non-responder      Pre_Responder 
#7524               2839               3203               2725 

#change the number of cells accordingly
feldman.pretreatment <- subset( feldman, 
                                cells=c( 
                                  base::sample( which( feldman@meta.data$treatment_response == "Pre_Non-responder" ), 2725, replace=FALSE ),
                                  base::sample( which( feldman@meta.data$treatment_response == "Pre_Responder" ), 2725, replace=FALSE ) ) )


DimPlot(feldman.pretreatment,split.by = "treatment_response", reduction= "umap", label=TRUE)

preprop_table <- table(feldman.pretreatment$seurat_clusters, feldman.pretreatment$treatment_response)
preprop_table
preprop_table[,1] <- preprop_table[,1]/sum(preprop_table[,1])
preprop_table[,2] <- preprop_table[,2]/sum(preprop_table[,2])
preprop_table
preprop <- as.data.frame((preprop_table))#.matrix((preprop_table))
which(preprop$Var1=="4")
preprop <- preprop[which(preprop$Var1=="4"),]
preprop

pdf("../feldman/feldman_all_therapyproportion.pdf")
ggplot(data=preprop, aes(x=Var1, y=Freq, fill = Var2)) +
  geom_bar(position="dodge",stat="identity", width=0.5) +
  ylab("Proportion") + 
  xlab("Cell type") +
  theme(text = element_text(size=20))+ guides(fill=guide_legend(title="Response"))
dev.off()
prop.test( c(246,190), c(2725,2725))

feldman.pretreatment.pre_No <- subset( feldman.pretreatment, 
                                       cells=which( feldman.pretreatment@meta.data$treatment_response == "Pre_Non-responder" ) )

feldman.pretreatment.pre_Yes <- subset( feldman.pretreatment, 
                                        cells=which( feldman.pretreatment@meta.data$treatment_response == "Pre_Responder" ) )


pdf("./feldman_UMAP.pdf")
DimPlot(feldman.pretreatment.pre_No, reduction = "umap")+ggtitle("Non-responder")
DimPlot(feldman.pretreatment.pre_Yes, reduction = "umap")+ggtitle("Responder")
DimPlot(feldman.pretreatment, reduction = "umap", split.by = "treatment_response")
dev.off()


#### Figure 7C: reanalysis of spatial data (zhang et al.)

# create seurat object

library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleCellExperiment)


dir.in <- '/fs/ess/PAS2579/NoJoon/GSE238264_HCC_spatial/GSE238264_RAW/'

SAMPLES <- c('HCC1R','HCC2R','HCC3R','HCC4R','HCC5NR','HCC6NR','HCC7NR')

for (i in 1:length(SAMPLES)){
  dir <- paste0(dir.in, SAMPLES[i])
  EXP <- Read10X(data.dir = paste0(dir,'/filtered_feature_bc_matrix'))
  IMAGE <-Read10X_Image(image.dir = paste0(dir,'/spatial'),slice = SAMPLES[i])
  coord <- IMAGE@coordinates
  coord <- coord[order(match(rownames(coord),colnames(EXP))),]
  identical(rownames(coord),colnames(EXP))
  DATA <- CreateSeuratObject(counts = EXP,project = SAMPLES[i],assay = 'Spatial')
  DATA@images$image <- new(Class = 'SlideSeq', assay = 'Spatial', key = 'image_', coordinates = coord[, 2:3])
  saveRDS(DATA, file = paste0(dir,'/',SAMPLES[i],'_seurat.rds'))
}



## non-responder: HCC6NR
i <-6

dir <- paste0(dir.in, SAMPLES[i])
dir
EXP <- Read10X(data.dir = paste0(dir,'/filtered_feature_bc_matrix'))
IMAGE <-Read10X_Image(image.dir = paste0(dir,'/spatial'),slice = SAMPLES[i])
coord <- IMAGE@coordinates
coord <- coord[order(match(rownames(coord),colnames(EXP))),]

CD8.COUNT <- EXP[which(rownames(EXP) %in% c('CD3E','CD3D','CD247','CD8B','CD8A','LCK')),]
CD8.COUNT.sum <- colSums(CD8.COUNT)
summary(CD8.COUNT.sum)

NK.COUNT <- EXP[which(rownames(EXP) %in% c('KLRD1','KLRF1','GNLY','CD7','NCR1','NKG7','KLRK1')),]

NK.COUNT.sum <- colSums(NK.COUNT)
summary(NK.COUNT.sum)

# CD8+ spots
IDX2 <- rep(0,length(CD8.COUNT.sum))
IDX2[which(CD8.COUNT.sum>=1)] <-1
IDX2 <- as.data.frame(IDX2)
colnames(IDX2) <- 'cut2'

IDX2$x <- coord[,2]
IDX2$y <- coord[,3]

ggplot(IDX2,aes(x,y,col= as.factor(cut2)))+geom_point(alpha=0.8)+theme_classic()+ scale_color_manual(values=c("grey", "red"))+
  ggtitle(paste0(SAMPLES[i],'_CD8','_cutoff=1'))+theme(legend.position="none")


# NK+ spots
IDX2 <- rep(0,length(NK.COUNT.sum))
IDX2[which(NK.COUNT.sum>=1)] <-1
IDX2 <- as.data.frame(IDX2)
colnames(IDX2) <- 'cut2'

IDX2$x <- coord[,2]
IDX2$y <- coord[,3]

ggplot(IDX2,aes(x,y,col= as.factor(cut2)))+geom_point(alpha=0.8)+theme_classic()+ scale_color_manual(values=c("grey", "blue"))+
  ggtitle(paste0(SAMPLES[i],'_NK','_cutoff=1'))+theme(legend.position="none")


## responder:HCC2R

i <- 2

dir <- paste0(dir.in, SAMPLES[i])
dir
EXP <- Read10X(data.dir = paste0(dir,'/filtered_feature_bc_matrix'))
IMAGE <-Read10X_Image(image.dir = paste0(dir,'/spatial'),slice = SAMPLES[i])
coord <- IMAGE@coordinates
coord <- coord[order(match(rownames(coord),colnames(EXP))),]

CD8.COUNT <- EXP[which(rownames(EXP) %in% c('CD3E','CD3D','CD247','CD8B','CD8A','LCK')),]
CD8.COUNT.sum <- colSums(CD8.COUNT)
summary(CD8.COUNT.sum)

NK.COUNT <- EXP[which(rownames(EXP) %in% c('KLRD1','KLRF1','GNLY','CD7','NCR1','NKG7','KLRK1')),]

NK.COUNT.sum <- colSums(NK.COUNT)
summary(NK.COUNT.sum)

IDX2 <- rep(0,length(CD8.COUNT.sum))
IDX2[which(CD8.COUNT.sum>=1)] <-1


IDX2 <- as.data.frame(IDX2)
colnames(IDX2) <- 'cut2'

IDX2$x <- coord[,2]
IDX2$y <- coord[,3]

# CD8+ spots
ggplot(IDX2,aes(x,y,col= as.factor(cut2)))+geom_point(alpha=0.8)+theme_classic()+ scale_color_manual(values=c("grey", "red"))+
  ggtitle(paste0(SAMPLES[i],'_CD8','_cutoff=1'))+theme(legend.position="none")


# NK+ spots
IDX2 <- rep(0,length(NK.COUNT.sum))
IDX2[which(NK.COUNT.sum>=1)] <-1


IDX2 <- as.data.frame(IDX2)
colnames(IDX2) <- 'cut2'

IDX2$x <- coord[,2]
IDX2$y <- coord[,3]

ggplot(IDX2,aes(x,y,col= as.factor(cut2)))+geom_point(alpha=0.8)+theme_classic()+ scale_color_manual(values=c("grey", "blue"))+
  ggtitle(paste0(SAMPLES[i],'_NK','_cutoff=1'))+theme(legend.position="none")


## ratio of CD8+ spots to NK+ spots across samples

RATIO <- vector()
for (i in 1:length(SAMPLES)){
  dir <- paste0(dir.in, SAMPLES[i])
  dir
  EXP <- Read10X(data.dir = paste0(dir,'/filtered_feature_bc_matrix'))
  IMAGE <-Read10X_Image(image.dir = paste0(dir,'/spatial'),slice = SAMPLES[i])
  coord <- IMAGE@coordinates
  coord <- coord[order(match(rownames(coord),colnames(EXP))),]

  CD8.COUNT <- EXP[which(rownames(EXP) %in% c('CD3E','CD3D','CD247','CD8B','CD8A','LCK')),]
  CD8.COUNT.sum <- colSums(CD8.COUNT)
  summary(CD8.COUNT.sum)

  NK.COUNT <- EXP[which(rownames(EXP) %in% c('KLRD1','KLRF1','GNLY','CD7','NCR1','NKG7','KLRK1')),]

  NK.COUNT.sum <- colSums(NK.COUNT)
  RATIO[i] <- length(which(CD8.COUNT.sum>=1))/length(which(NK.COUNT.sum>=1))
  
}

DF <- data.frame(ratio = RATIO)
DF$group <- rep(c('R','NR'), c(4,3))
ggplot(DF, aes(x = group, y = ratio)) +  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1)) + theme_classic()


## Figure 7D: reanalysis of scRNA-seq data (zheng et al., pancancer)

CD8_integrated <- readRDS('/fs/ess/PAS2579/NoJoon/pancancer/int.CD8.S35.sce.merged.rds')
UMAP.harmony <- as.data.frame(reducedDims(CD8_integrated)[['harmony.umap']])
UMAP.harmony$group <- as.character(CD8_integrated$meta.cluster)

ggplot(UMAP.harmony, aes(harmony.umap_1, harmony.umap_2, color= group))+ geom_point() + theme_classic()

GENE <- read.csv('/fs/ess/PAS2579/NoJoon/geneList2_song.csv')  # alternative set2
IL2 <- read.csv('/fs/ess/PAS2579/NoJoon/pancancer/IL2_pathway.csv')
IFNalpha <- read.csv('/fs/ess/PAS2579/NoJoon/pancancer/IFNalpha_pathway.csv')
RAR.target <- read.csv('/fs/ess/PAS2579/NoJoon/set3.csv',header=T)    # IRIS3 RAR target


IL2 <- IL2$Symbol
IL2 <- IL2[nzchar(IL2)]

IFNalpha <- IFNalpha$Symbol
IFNalpha <- IFNalpha[nzchar(IFNalpha)]

## if IRIS3 RAR
RAR <- RAR.target$Set3_RAR


count <- CD8_integrated@assays$data$exprs

length(intersect(IL2,toupper(GENE$gene)))  # 0
length(intersect(IFNalpha,toupper(GENE$gene)))  # 0
length(intersect(toupper(RAR),toupper(GENE$gene)))  # 0


TARGET.GENE <- c(GENE$gene,IL2,IFNalpha,toupper(RAR))

TARGET.count <- count[which(toupper(rownames(count))%in% TARGET.GENE),]
TARGET.count <- t(TARGET.count)



a <- vector()

IND <- which(colnames(TARGET.count)%in% toupper(GENE$gene))

for(i in 1:ncol(TARGET.count)){
  if (i %in% IND){
    temp <- colnames(TARGET.count)[i]
    tt <- GENE[which(toupper(GENE$gene)==temp),]$exp
    a[i] <- ifelse(tt=='High',0.1,0.1)
  }else{
    a[i] <- 0.1
  }
}



a <- as.matrix(a)
TARGET.score <- TARGET.count %*% a




summary(TARGET.score[,1])  
hist(TARGET.score[,1])

TT <-rownames(TARGET.count)[-which(is.na(TARGET.score))]

TARGET.score <- TARGET.score[-which(is.na(TARGET.score))]
names(TARGET.score) <- TT
UMAP.sub <- UMAP.harmony[which(rownames(UMAP.harmony)%in% names(TARGET.score)),]


TARGET.score.mod <- TARGET.score

for (i in 1:length(TARGET.score.mod)){
  TARGET.score.mod[i] <- ifelse(TARGET.score.mod[i]>2,2,TARGET.score.mod[i])
  TARGET.score.mod[i] <- ifelse(TARGET.score.mod[i]<(-2),(-2),TARGET.score.mod[i])
}

hist(TARGET.score.mod)


df <- data.frame(UMAP.sub,TARGET.score,TARGET.score.mod)

mid <- median(TARGET.score.mod)
ggplot(df,aes( harmony.umap_1, harmony.umap_2,color=TARGET.score.mod))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()


## extract the signature from pancancer data
temp <- CD8_integrated@assays$data$exprs
CV <- apply(temp,1,sd)
INDS <- which(is.na(CV))

CD8_sub <- CD8_integrated[-INDS,]
stats <- modelGeneCV2(CD8_sub,assay.type = 'exprs')
HIGH <- getTopHVGs(stats,var.field = 'ratio')
FINAL <- union(HIGH,TARGET.GENE)
FINAL <- FINAL[which(FINAL%in% rownames(CD8_sub))]
CD8_sub <- CD8_sub[FINAL,]

left <- setdiff(1:nrow(df),INDEX)
library(scran)

RST <- list()
for (i in 10001:10020){
  set.seed(i)
  left.sub <- sample(left,1000)
  DF.sub <- DF[c(left.sub,INDEX),]
  #CD8.sub <- CD8_integrated[,which(colnames(CD8_integrated)%in% rownames(DF.sub))]
  CD8.sub <- CD8_sub[,which(colnames(CD8_sub)%in% rownames(DF.sub))]
  DF.sub <-DF.sub[order(match(rownames(DF.sub),colnames(CD8.sub))),]
  
  #LOW <- findMarkers(CD8.sub,groups = DF.sub$A,assay.type = 'exprs')
  #LOW <- findMarkers(CD8.sub,groups = as.factor(DF.sub$A),assay.type = 'exprs',test.type='t')
  LOW <-scoreMarkers(CD8.sub,groups = DF.sub$A,assay.type = 'exprs')
  LOW[[1]][order(LOW[[1]]$mean.AUC, decreasing=TRUE),1:4]   # this is the upregulated genes in group1
  RST[[i-10000]] <- LOW
  
}

saveRDS(RST,file = '/fs/ess/PAS2579/NoJoon/pancancer/low.rds')


low <-readRDS('/fs/ess/PAS2579/NoJoon/pancancer/low.rds')


MARKERS <- function(x){
  temp <- x[[1]]
  t2 <-temp[order(temp$mean.AUC, decreasing=TRUE),] 
  rownames(t2[1:50,])
}

low.marker <-lapply(low,MARKERS)
markers1 <- Reduce(intersect, low.marker)


c10 <- readRDS('/fs/ess/PAS2579/NoJoon/pancancer/high_C10.rds')

c10.marker <-lapply(c10,MARKERS)
markers.c10 <- Reduce(intersect, c10.marker)

c789 <- readRDS('/fs/ess/PAS2579/NoJoon/pancancer/high_C7-9.rds')

c789.marker <-lapply(c789,MARKERS)
markers.c789 <- Reduce(intersect, c789.marker)


####### pre
yost.pre.CD8 <- readRDS('/fs/ess/PAS2579/NoJoon/Yost/yost.pre.CD8.rds')
responder <- c('su001','su002','su003','su004','su009','su012')
noResponder <- c('su005','su006','su007','su008','su010')


yesCell <- colnames(yost.pre.CD8)[which(yost.pre.CD8$patient %in% responder)] 
NoCell <- colnames(yost.pre.CD8)[which(yost.pre.CD8$patient %in% noResponder)] 

count.data <- GetAssayData(yost.pre.CD8)

## low
low.count <- count.data[which(rownames(count.data)%in% markers1),]
low.count <- t(low.count)
a <- rep(1,ncol(low.count))
a <- as.matrix(a)
low.score <- low.count%*%a

summary(low.score[,1])  # there are ~700 NaN
hist(low.score[,1])

low.score.mod <- low.score[,1]
hist(low.score.mod)

for (i in 1:length(low.score.mod)){
  #low.score.mod[i] <- ifelse(low.score.mod[i]>45,45,low.score.mod[i])
  low.score.mod[i] <- ifelse(low.score.mod[i]<20,20,low.score.mod[i])
}

hist(low.score.mod)




UMAP <- yost.pre.CD8@reductions$umap@cell.embeddings

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(low.score)),]

df <- data.frame(UMAP.sub,low.score,low.score.mod)
RESPONSE <- rep(-1,nrow(df))
for (i in 1:nrow(df)){
  RESPONSE[i] <- ifelse(rownames(df)[i] %in% yesCell,1,0)
}

df$responder <- RESPONSE


ggplot(df, aes(x=as.factor(RESPONSE), y=low.score.mod)) + 
  geom_boxplot()+theme_classic()+ggtitle('Yost_pre_low')


##
c10.count <- count.data[which(rownames(count.data)%in% markers.c10),]
c10.count <- t(c10.count)
a <- rep(1,ncol(c10.count))
a <- as.matrix(a)
c10.score <- c10.count%*%a

summary(c10.score[,1])  # there are ~700 NaN
hist(c10.score[,1])

c10.score.mod <- c10.score[,1]
hist(c10.score.mod)

for (i in 1:length(c10.score.mod)){
  c10.score.mod[i] <- ifelse(c10.score.mod[i]>80,80,c10.score.mod[i])
  c10.score.mod[i] <- ifelse(c10.score.mod[i]<20,20,c10.score.mod[i])
}

hist(c10.score.mod)




UMAP <- yost.pre.CD8@reductions$umap@cell.embeddings

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(c10.score)),]

df <- data.frame(UMAP.sub,c10.score,c10.score.mod)
RESPONSE <- rep(-1,nrow(df))
for (i in 1:nrow(df)){
  RESPONSE[i] <- ifelse(rownames(df)[i] %in% yesCell,1,0)
}

df$responder <- RESPONSE

ggplot(df, aes(x=as.factor(RESPONSE), y=c10.score.mod)) + 
  geom_boxplot()+theme_classic()+ggtitle('Yost_pre_c10')

df.sub <- df[which(df$c10.score>=37&df$c10.score<=65),]

ggplot(df.sub, aes(x=as.factor(responder), y=c10.score.mod)) + 
  geom_boxplot()+theme_classic()+ggtitle('Yost_pre_c10')

t.test(df.sub[which(df.sub$responder==1),]$c10.score.mod,df.sub[which(df.sub$responder==0),]$c10.score.mod,var.equal = T)

#############################################################################################
############## Figure S11C: scRNA-seq DATA (Daniel et al)

library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)

DAT <- Read10X(data.dir = '/fs/ess/PAS2579/NoJoon/GSE188666/out',gene.column=1)
DAT <- CreateSeuratObject(counts = DAT,project = "GSE18866",min.cells = 0, min.features = 0)

a <-sapply(strsplit(colnames(DAT),'_',1),'[',3)
table(a)  # contains D8 and D21

DAT$time <- a


meta <- read.table('/fs/ess/PAS2579/NoJoon/GSE188666/out/GSE188666_scRNA_LCMV_metadata.tsv')
table(meta$tissue)   # liver, spleen and lung

identical(colnames(DAT),rownames(meta))  # check if the order of cells are the same

DAT <- AddMetaData(DAT, meta)

Idents(DAT) <- DAT$ident
UMAP <- data.frame(UMAP_1=DAT$UMAP_1,UMAP_2=DAT$UMAP_2)
UMAP <- as.matrix(UMAP)
DAT[['umap']] <-CreateDimReducObject(embeddings = UMAP,key = 'UMAP_',assay = DefaultAssay(DAT))
DAT$clusters <- DAT$ident


DimPlot(DAT,label=T)



count <- DAT@assays$RNA@counts
UMAP <- DAT@reductions$umap@cell.embeddings

### NKRT signature
GENE <- read.csv('/fs/ess/PAS2579/NoJoon/geneList2_song.csv')   # alternative set2
GENE.count <- count[which(toupper(rownames(count))%in% toupper(GENE$gene)),]
GENE.count <- t(GENE.count)

a <- rep(1,ncol(GENE.count))
a <- as.matrix(a)
GENE.score <- GENE.count%*%a

summary(GENE.score[,1])  
hist(GENE.score[,1])
hist(log(GENE.score[,1]+1))


rownames(GENE.score) <- rownames(GENE.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(GENE.score)),]

df <- data.frame(UMAP.sub,GENE.score, log(GENE.score+1))

mid <- median(log(GENE.score))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(GENE.score)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()



## RAR

RAR.target <- read.csv('/fs/ess/PAS2579/NoJoon/set3.csv',header=T)

RAR.count <- count[which(rownames(count)%in% RAR.target$Set3_RAR),]
RAR.count <- t(RAR.count)

a <- rep(1,ncol(RAR.count))
a <- as.matrix(a)
RAR.score <- RAR.count%*%a

summary(RAR.score[,1])


rownames(RAR.score) <- rownames(RAR.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(RAR.score)),]

df <- data.frame(UMAP.sub,RAR.score)
mid <- median(log(RAR.score))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(RAR.score)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()


## IL2 
IL2 <- read.csv('/fs/ess/PAS2579/NoJoon/pancancer/IL2_pathway.csv')

IL2 <- IL2$Symbol
IL2 <- IL2[nzchar(IL2)]

IL2.count <- count[which(toupper(rownames(count))%in% IL2),]
IL2.count <- t(IL2.count)

a <- rep(1,ncol(IL2.count))
a <- as.matrix(a)
IL2.score <- IL2.count%*%a

summary(IL2.score[,1])


rownames(IL2.score) <- rownames(IL2.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(IL2.score)),]

df <- data.frame(UMAP.sub,IL2.score)
mid <- median(log(IL2.score+1))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(IL2.score+1)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()





## IFNalpha
IFNalpha <- read.csv('/fs/ess/PAS2579/NoJoon/pancancer/IFNalpha_pathway.csv')
IFNalpha <- IFNalpha$Symbol
IFNalpha <- IFNalpha[nzchar(IFNalpha)]

IFN.count <- count[which(toupper(rownames(count))%in% IFNalpha),]
IFN.count <- t(IFN.count)

a <- rep(1,ncol(IFN.count))
a <- as.matrix(a)
IFN.score <- IFN.count%*%a

summary(IFN.score[,1])


rownames(IFN.score) <- rownames(IFN.count)

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(IFN.score)),]

df <- data.frame(UMAP.sub,IFN.score)
mid <- median(log(IFN.score+1))
ggplot(df,aes(UMAP_1,UMAP_2,color=log(IFN.score+1)))+geom_point(size=0.5) +
  scale_colour_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()


##### combine RAR and NKRT together


GENE <- read.csv('/fs/ess/PAS2579/NoJoon/geneList2_song.csv')  # alternative set2

TARGET.GENE <- c(GENE$gene,IL2,IFNalpha,toupper(RAR))

TARGET.count <- count[which(toupper(rownames(count))%in% TARGET.GENE),]
TARGET.count <- t(TARGET.count)


a <- vector()

IND <- which(colnames(TARGET.count)%in% toupper(GENE$gene))

for(i in 1:ncol(TARGET.count)){
  if (i %in% IND){
    temp <- colnames(TARGET.count)[i]
    tt <- GENE[which(toupper(GENE$gene)==temp),]$exp
    a[i] <- ifelse(tt=='High',0.1,0.1)
  }else{
    a[i] <- 0.1
  }
}



a <- as.matrix(a)
TARGET.score <- TARGET.count %*% a

UMAP.sub <- UMAP[which(rownames(UMAP)%in% rownames(TARGET.score)),]
TARGET.score.mod <- log(TARGET.score[,1])


df <- data.frame(UMAP.sub,TARGET.score,TARGET.score.mod)

mid <- median(TARGET.score.mod)
ggplot(df,aes(UMAP_1, UMAP_2,color=TARGET.score.mod))+geom_point(size=0.2) +
  scale_colour_gradient2(low = "royalblue",
                         mid = "white",
                         high = "red",
                         midpoint = mid)+theme_classic()


#############################################################################################
#Griffths - Supp figure S13

library(tidyverse)
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(readr)
library(data.table)
library(cowplot)
library(readxl)
#library(biomaRt)

# data directory (change on your computer)
# For PBMC Seurat tutorial data, put "filtered_gene_bc_matrices/" within the data_dir
data_dir <- c("./")
# PBMC Seurat tutorial data

grif.data <- read.table(file = paste0(data_dir,"GSE130157.cell_annotations.txt"), sep="\t", header=TRUE, row.names=1)
grif.data <- grif.data %>% 
  mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
  na.omit() ##removing blanks

grif.data$Time.Point <- gsub("C1", "Pre", grif.data$Time.Point)
grif.data$Time.Point <- gsub("C.*", "Post", grif.data$Time.Point)


grif.data1 <- read.table(file = paste0(data_dir,"14546.14554.14562.RawCounts.txt"), sep="\t", header=TRUE, row.names = 1)
grif.data2 <- read.table(file = paste0(data_dir,"15424R.15435R.RawCounts.txt"), sep="\t", header=TRUE, row.names = 1)
grif.data1 <- grif.data1[3:52271]
grif.data2 <- grif.data2[3:18517]

for (i in 1:52269){
  colnames(grif.data1)[i] <-sub('.','', colnames(grif.data1)[i])
}
for (i in 1:18515){
  colnames(grif.data2)[i] <- sub('.','', colnames(grif.data2)[i])
}

header1 <- rownames(grif.data)[rownames(grif.data) %in% colnames(grif.data1)]
header2 <- rownames(grif.data)[rownames(grif.data) %in% colnames(grif.data2)]

grif.data1 <- grif.data1 %>% dplyr::select(header1)
grif.data2 <- grif.data2 %>% dplyr::select(header2)

grif.data0 <- cbind(grif.data1, grif.data2)
#save(grif.data, grif.data0, file = "./grifcountmeta.RData")
load("./grifcountmeta.RData")

#preprocessing for removing NA in subcluster


#ensemble to official gene
ensembl_ids <- rownames(grif.data0)
# 3. Get a list of 'official gene symbol'
ens <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes.tbl <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                   filters="ensembl_gene_id",
                   values=ensembl_ids, mart=ens)
genes.tbl <- genes.tbl[!duplicated(genes.tbl$ensembl_gene_id),]
# 4. Convert the ensembl to symbol
#grif.dataT0 <- t(grif.data0)
#header0 <- colnames(grif.dataT0)[colnames(grif.dataT0) %in% genes.tbl$ensembl_gene_id]
#grif.dataT0 <- grif.dataT0[,header0]
#grif.dataT0$X <- colnames(n) 

grif.data0$X <- rownames(grif.data0)
merge.tbl <- merge(x=grif.data0, y=genes.tbl,
                   by.x="X", by.y="ensembl_gene_id",
                   all.x=TRUE, all.y=FALSE)

merge.tbl <- merge.tbl[,c(1,56898,seq(2,56897))]
merge.tbl <- merge.tbl[-which(is.na(merge.tbl$hgnc_symbol)),]
merge.tbl <- merge.tbl[!duplicated(merge.tbl[,2]),]
merge.tbl <- merge.tbl[-c(1588),] #empty feature
rownames(merge.tbl) <- merge.tbl[,2]
merge.tbl <- merge.tbl[,c(3:56898)]
# for (i in 1:22953){
#   rownames(merge.tbl)[i] <- gsub('-','', rownames(merge.tbl)[i])
# }
write.csv(colnames(merge.tbl), "./colnames_griffths.csv", row.names = TRUE)
write.csv(rownames(merge.tbl), "./rownames_griffths.csv", row.names = TRUE)

# Initialize the Seurat object with the raw (non-normalized data).
grif <- CreateSeuratObject(counts = grif.data0, 
                           project = "NKcells_count", 
                           meta.data = grif.data,
                           min.cells = 3, #low quality genes
                           min.features = 200) #at least 200 features

save(grif, file=paste0(data_dir,"./grif.RData"))
save(grif, file = paste0(data_dir,"./grif_convert.RData"))
load(paste0(data_dir,"./grif_convert.RData"))

grif$treatment_response <- paste(grif@meta.data$Time.Point, grif@meta.data$Responder, sep="_")

# Normalize counts
grif <- NormalizeData(grif)

# Identify the most informative genes
grif <- FindVariableFeatures(grif, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# Scale the genes
grif <- ScaleData(grif)
# Compute PCA (50)
grif <- RunPCA(grif, features = VariableFeatures(object = grif))
#rectangular data

#From here graph data
# Construct cell-cell nearest neighbors graph
grif <- FindNeighbors(grif, dims = 1:10)
# Cluster the graph using Louvain
grif <- FindClusters(grif, resolution = 0.5)

# # save counts for singleR
counts <- grif[["RNA"]]@counts
metadata <- grif@meta.data
save(counts,metadata,file = "./grifcounts_norm.RData")
save(grif, file="./grifafternorm.RData")
save(counts,metadata, file=paste0(data_dir, "./grifcounts_norm_conver.RData"))
save(grif, file= paste0(data_dir,"./grifafternorm_convert.RData"))
save(counts,metadata, file=paste0(data_dir, "./grifcounts_norm_subcluster.RData"))
save(grif, file= paste0(data_dir,"./grifafternorm_subcluster.RData"))
save(counts,metadata, file=paste0(data_dir, "./grifcounts_norm_conver_subcluster.RData"))
save(grif, file= paste0(data_dir,"./grifafternorm_convert_subcluster.RData"))
#load("./grifafternorm.RData")
#load(paste0(data_dir,"./grifafternorm_convert.RData"))
load("../../griffiths/grifafternorm_convert_subcluster.RData")
#load("./griffiths/grifafternorm_convert_subcluster.RData")
# Compute UMAP coordinates for visualization
grif <- RunUMAP(grif, dims = 1:10) #only visualization usage. UMAP  nonlinear / using 10 pca features /hyper parameter


#Reference based annotation

load('./griffiths/grifcounts_norm_conver_subcluster.RData')
sce <- SingleCellExperiment(assays = list(counts = counts))

# lognorm transform
sce <- logNormCounts(sce)

# load various reference databases
ref_hum <- HumanPrimaryCellAtlasData()
ref_dat <- DatabaseImmuneCellExpressionData()
ref_nov <- NovershternHematopoieticData()
ref_mon <- MonacoImmuneData()

# cell type prediction
pred_hum <- SingleR(test = sce, ref = ref_hum, labels = ref_hum$label.main)
pred_dat <- SingleR(test = sce, ref = ref_dat, labels = ref_dat$label.main)
pred_nov <- SingleR(test = sce, ref = ref_nov, labels = ref_nov$label.main)
pred_mon <- SingleR(test = sce, ref = ref_mon, labels = ref_mon$label.main)

pred_list <- list(pred_hum,
                  pred_dat,
                  pred_nov,
                  pred_mon)
save(pred_list,
     file = "./grif_norm_singleR_predictions_subcluster.RData")


#load("./grif_norm_singleR_predictions.RData")

#prediction diagnostics - HumanPrimaryCellAtlasData
pred_plot_hum <- plotScoreHeatmap(pred_list[[1]],
                                  annotation_col=as.data.frame(metadata[,"Sub.Cluster",
                                                                        drop=FALSE]))
ggsave(plot = pred_plot_hum, 
       filename = "./grif_norm_celltype_pred_heatmap_human_subcluster.pdf", 
       device = "pdf",
       height  = 8,
       width  = 9)

pred_table_hum <- table(pred_list[[1]]$labels, metadata[,"Sub.Cluster"])
write.csv(pred_table_hum, 
          file = "./grif_norm_celltype_pred_table_human_subcluster.csv",
          quote = FALSE,
          row.names = TRUE)


# prediction diagnostics - DatabaseImmuneCellExpressionData
pred_plot_dat <- plotScoreHeatmap(pred_list[[2]], 
                                  annotation_col=as.data.frame(metadata[,"Sub.Cluster",
                                                                        drop=FALSE]))
ggsave(plot = pred_plot_dat, 
       filename = "./grif_norm_celltype_pred_heatmap_database_subcluster.pdf", 
       device = "pdf",
       height  = 8,
       width  = 9)

pred_table_dat <- table(pred_list[[2]]$labels, metadata[,"Sub.Cluster"])
write.csv(pred_table_dat, 
          file = "./grif_norm_celltype_pred_table_database_subcluster.csv",
          quote = FALSE,
          row.names = TRUE)


#prediction diagnostics - NovershternHematopoieticData
pred_plot_nov <- plotScoreHeatmap(pred_list[[3]], 
                                  annotation_col=as.data.frame(metadata[,"Sub.Cluster",
                                                                        drop=FALSE]))
ggsave(plot = pred_plot_nov, 
       filename = "./grif_norm_celltype_pred_heatmap_novershtern_subcluster.pdf", 
       device = "pdf",
       height  = 8,
       width  = 9)

pred_table_nov <- table(pred_list[[3]]$labels, metadata[,"Sub.Cluster"])
write.csv(pred_table_nov, 
          file = "./grif_norm_celltype_pred_table_novershtern_subcluster.csv",
          quote = FALSE,
          row.names = TRUE)


# prediction diagnostics - MonacoImmuneData
pred_plot_mon <- plotScoreHeatmap(pred_list[[4]],
                                  annotation_col=as.data.frame(metadata[,"Sub.Cluster",
                                                                        drop=FALSE]))
ggsave(plot = pred_plot_mon, 
       filename = "./grif_norm_celltype_pred_heatmap_monaco_subcluster.pdf", 
       device = "pdf",
       height  = 8,
       width  = 9)

pred_table_mon <- table(pred_list[[4]]$labels, metadata[,"Sub.Cluster"])
write.csv(pred_table_mon, 
          file = "./grif_norm_celltype_pred_table_monaco_subcluster.csv",
          quote = FALSE,
          row.names = TRUE)


table(grif@meta.data$treatment_response)
#Post_Non.Responder     Post_Responder  Pre_Non.Responder      Pre_Responder 
#17310                  19276           8848                   11462 

grif.pretreatment <- subset( grif, 
                             cells=c( 
                               base::sample( which( grif@meta.data$treatment_response == "Pre_Non.Responder" ), 8848, replace=FALSE ),
                               base::sample( which( grif@meta.data$treatment_response == "Pre_Responder" ), 8848, replace=FALSE ) ) )

preprop_table <- table(grif.pretreatment$Sub.Cluster, grif.pretreatment$treatment_response)
preprop_table
preprop_table[,1] <- preprop_table[,1]/sum(preprop_table[,1])
preprop_table[,2] <- preprop_table[,2]/sum(preprop_table[,2])
preprop_table
preprop <- as.data.frame((preprop_table))#.matrix((preprop_table))
preprop_act <- preprop[which(preprop$Var1=="NK_Cell_Activated"),]
preprop_rest <- preprop[which(preprop$Var1=="NK_Cell_Resting"),]

# Supp - Fig 13B

pdf("./grif_proportion.pdf")
ggplot(data=preprop_act, aes(x=Var1, y=Freq, fill = Var2)) +
  geom_bar(position="dodge",stat="identity", width=0.5) + ggtitle("NK cell activated")+
  ylab("Proportion") + 
  xlab("Cell type") +
  theme(text = element_text(size=20))+ guides(fill=guide_legend(title="Response"))

ggplot(data=preprop_rest, aes(x=Var1, y=Freq, fill = Var2)) +
  geom_bar(position="dodge",stat="identity", width=0.5) + ggtitle("NK cell resting")+
  ylab("Proportion") + 
  xlab("Cell type") +
  theme(text = element_text(size=20))+ guides(fill=guide_legend(title="Response"))

dev.off()
#activated
prop.test( c(198,127), c(8848,8848))
#resting
prop.test( c(1479,805), c(8848,8848))


grif.pretreatment.pre_No <- subset( grif.pretreatment, 
                                    cells=which( grif.pretreatment@meta.data$treatment_response == "Pre_Non.Responder" ) )

grif.pretreatment.pre_Yes <- subset( grif.pretreatment, 
                                     cells=which( grif.pretreatment@meta.data$treatment_response == "Pre_Responder" ) )

Idents(grif.pretreatment) <- "Sub.Cluster"
Idents(grif.pretreatment.pre_No) <- "Sub.Cluster"
Idents(grif.pretreatment.pre_Yes) <- "Sub.Cluster"


#Supp - Fig 13A
pdf("./grif_UMAP_nolegend.pdf")
DimPlot(grif.pretreatment.pre_No, reduction = "umap")+NoLegend()+ggtitle("Non-responder")
DimPlot(grif.pretreatment.pre_Yes, reduction = "umap")+NoLegend()+ggtitle("Responder")
DimPlot(grif.pretreatment, reduction = "umap", split.by = "treatment_response")+NoLegend()
dev.off()





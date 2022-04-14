# "Single-cell gene expression analysis of cryopreserved equine bronchoalveolar cells"
# Sophie E. Sage, Pamela Nicholson, Laureen M. Peters, Tosso Leeb, Vidhya Jagannathan*, Vinzenz Gerber*,
# *These authors contributed equally to this work.

# Set up working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set up R environment
library(Seurat)
library(clustree)
library(AnnotationHub)
library(dplyr)
library(Matrix)
library(gdata)
library(biomaRt)
library(stringr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(reshape2)
library(data.table)
library(dittoSeq)
library(ggpubr)

# PRE-PROCESSING
# Raw sequencing data (fastq files) was converted to a count matrix of gene expression values using the Cell Ranger (v6.0) standard workflow. 
# The 3'-untranslated regions of the genes in the reference genome (Equus caballus NCBI annotation release 103) were manually extended by 2 kb.  
# Gene IDs were converted to gene symbols using AnnotationHub() before running CellRanger
# 'sc2021_ncbi_extend.RData' RData contains the output of CellRanger 6.0 

load('sc2021_ncbi_extend.RData') 
sc21 <- sc2021 # "sc2021" is the raw data Seurat object


# QUALITY CONTROL & FILTERING

# Add mitochondrial reads to the Seurat object's meta.data
mito.genes<-c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB","MT-tRNA-Phe",  
              "MT-s-rRNA","MT-tRNA-Val","MT-l-rRNA","MT-tRNA-Leu","MT-tRNA-Ile","MT-tRNA-Gln","MT-tRNA-Met","MT-tRNA-Trp","MT-tRNA-Ala",
              "MT-tRNA-Cys","MT-tRNA-Tyr","MT-tRNA-Ser","MT-tRNA-Asp","MT-tRNA-Lys","MT-tRNA-Gly","MT-tRNA-Arg","MT-tRNA-His","MT-tRNA-Ser.1",
              "MT-tRNA-Leu.1","MT-tRNA-Glu","MT-tRNA-Thr","MT-tRNA-Pro")#36/37 MT genes because tRNA-Asn is absent from this dataset
sc21[["percent.mito"]] <- PercentageFeatureSet(sc21, features=mito.genes)

# Relationship between number of UMIs and percentage of mitochondrial genes
FeatureScatter(sc21, feature1 = "nCount_RNA", feature2 = "percent.mito") 

# Data visual inspection
table(sc21$orig.ident) #number of cells per sample => 5,408 cells
VlnPlot(sc21, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
VlnPlot(sc21, features = "nFeature_RNA")
VlnPlot(sc21, features = "percent.mito")
FeatureScatter(sc21, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(sc21, feature1 = "nCount_RNA", feature2 = "percent.mito")

# Supplementary Figure 1: Cell filtering based on feature count and percentage of mitochondrial reads
Plot_FilterVln <- ((VlnPlot(sc21, features = "nFeature_RNA") + 
                      geom_hline(yintercept=200, linetype="dashed", color = "red") +
                      geom_hline(yintercept=6500, linetype="dashed", color = "red") +
                      labs(title= NULL, x=NULL, y="Feature count"))+
                     (VlnPlot(sc21, features = "percent.mito") + 
                        geom_hline(yintercept=15, linetype="dashed", color = "red") +
                        labs(title= NULL, x=NULL, y="Mitochondrial reads %")) +
                     plot_layout(ncol = 1)) &
  (theme(axis.title = element_text(size=10), axis.text.x = element_blank(), 
         axis.ticks.x = element_blank(), axis.text.y = element_text(size=8))) &
  scale_fill_discrete(limits = c("ss1", "ss2", "ss3"),
                      labels = c("Horse 1", "Horse 2", "Horse 3"))

Plot_FilterFeat <- ((FeatureScatter(sc21, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
                       geom_hline(yintercept=200, linetype="dashed", color = "red") +
                       geom_hline(yintercept=6500, linetype="dashed", color = "red") + 
                       labs(title= NULL, x="RNA count", y="Feature count")) +
                      (FeatureScatter(sc21, feature1 = "nCount_RNA", feature2 = "percent.mito") + 
                         geom_hline(yintercept=15, linetype="dashed", color = "red") 
                       + labs(title= NULL, x="RNA count", y="Mitochondrial reads %")) +
                      plot_layout(ncol = 1)) &
  theme(axis.title = element_text(size=10), axis.text = element_text(size=8)) &
  NoLegend()

dpi=300
tiff(file='sc21_Filters.tiff', width = dpi*9, height = dpi*6, units = "px",res = dpi)
(Plot_FilterFeat | Plot_FilterVln) + plot_annotation(tag_levels = 'A')
dev.off()

# Based on the scatter and violin plots, the following filtering thresholds are selected:
sc21_flt <- subset(sc21, nFeature_RNA >= 200 & nFeature_RNA<6500) #5,256 cells (152 cells filtered)
sc21_flt <- subset(sc21_flt, percent.mito <= 15) #4,631 cells (625+152=777=14% filtered)
sc21 <- sc21_flt


# NORMALIZATION

sc21 <- NormalizeData(sc21, normalization.method = "LogNormalize", scale.factor = 10000)


#VARIABLE FEATURES SELECTION

# Variable features selection: we calculate a subset of features that exhibit high cell-to-cell variation in the dataset
sc21 <- FindVariableFeatures(sc21, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(sc21))
vf_top10 <- head(VariableFeatures(sc21), 10)
vf_plot <- VariableFeaturePlot(sc21)
LabelPoints(plot = vf_plot,
            points = vf_top10, repel = TRUE, xnudge = 0, ynudge = 0)


# SCALING

sc21 <- ScaleData(object = sc21)


# QUALITY CONTROL

# Cell cycle analysis
# Convert equine genes to human orthologues 
convertHumanGeneList <- function(x){
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest")
  horse = useEnsembl("ensembl", dataset = "ecaballus_gene_ensembl", mirror = "uswest")
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("external_gene_name"), martL = horse, uniqueRows=T)
  humanx <- unique(genes[, 2])
  print(head(humanx))
  return(humanx)
}

# cc.genes.updated.2019 is a list of genes used in cell-cycle regression, updated with 2019 symbols
# It is a list of two vectors: s.genes (Genes associated with S-phase) and g2m.genes (Genes associated with G2M-phase)
eq.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
eq.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
sc21 <- CellCycleScoring(sc21, s.features = eq.s.genes, g2m.features = eq.g2m.genes, set.ident = TRUE)
table(Idents(sc21)) # the previous step fixed the cell cycle phase as default identities in Seurat
sc21 <- SetIdent(sc21, value = sc21$orig.ident) # change Id back to sample origin


# DIMENSIONALITY REDUCTION (PCA)

sc21 <- RunPCA(sc21)
DimPlot(sc21, reduction = "pca")
DimPlot(sc21, reduction = "pca", group.by = "Phase") 
# Cells are not clustering based on cell cycle phase; thus regression is not necessary.
head(Loadings(sc21, reduction = "pca")[, 1:10])
# Selection of the number of PCs...
#...first visually with an elbow plot
ElbowPlot(sc21, ndims = 40) +
  geom_vline(xintercept = 13, linetype="dashed", color = "blue") +
  geom_vline(xintercept = 16, linetype="dashed", color = "green") #16 PCs look appropriate
#...and with a quantitative technique
pct <- sc21[["pca"]]@stdev / sum(sc21[["pca"]]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu>90 & pct<5)[1]
co1 #42
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 #10
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct, 
                           cumu = cumu, 
                           rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw() + ggtitle("Elbow plot (quantitative)")
# We select 16 PCs because UMAP visualization better differentiates cell clusters when using 16 PCs compared to 10 PCs (not shown).
DimHeatmap(sc21, dims = 1:16, cells = 500, balanced = TRUE) 
VizDimLoadings(sc21, dims = 16, ncol = 1) + theme_minimal(base_size = 8) 


# CLUSTERS VISUALIZATION WITH UMAP

sc21 <- RunUMAP(sc21, dims = 1:16) 
DimPlot(sc21, reduction = "umap")
DimPlot(sc21, reduction = "umap", group.by = "Phase")


# CELL CLUSTERING

sc21 <- FindNeighbors(sc21, dims = 1:16) 
sc21 <- FindClusters(sc21, resolution = seq(0.1, 1.2, by=0.1))
head(sc21@meta.data)

# Clustering resolution is chosen based on 2 criteria:
# 1) Cluster stability (assessed with the clustree package)
# 2) Visual fit to the data set (assessed with UMAP)
# Visualize how clusters sub-divide at increasing resolution:
clustree(sc21@meta.data[,grep("RNA_snn_res", colnames(sc21@meta.data))],
                   prefix = "RNA_snn_res.")
# Visualize the UMAP at different resolutions:
(DimPlot(object = sc21, group.by=grep("res",colnames(sc21@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

(DimPlot(object = sc21, group.by=grep("res",colnames(sc21@meta.data),value = TRUE)[9:12], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))

DimPlot(sc21, reduction = "umap", group.by = "RNA_snn_res.1", label=T)
# We choose resolution 1.0 (16 clusters)


#CELL ANNOTATION

sc21 <- SetIdent(sc21, value = sc21$RNA_snn_res.1) # set up resolution = 1.0

# Create vectors of marker genes(features) for cell types 
feat_mac <- c("CD68","CD163") # markers for monocytes-macrophages (Mo/Ma)
feat_T <- c("CD2", "CD3D","CD3E","CD3G") # markers for T cells
feat_B <- c("MS4A1", "CD79A","CD79B") # markers for B/Plasma cells
feat_neut <- c("TG","RGS2","ILT11B","CSF3R") # markers for neutrophils
feat_mast <- c("LTC4S","HPGDS","GCSAML","MS4A2") # markers for mast cells
feat_DC <- c("CD83","CCR7")# markers for dendritic cells 


# Plot the expression of cell type specific features
FeaturePlot(sc21, features = feat_mac, ncol=2)
VlnPlot(sc21, features = feat_mac) 

FeaturePlot(sc21, features = feat_T, ncol=2)
VlnPlot(sc21, features = feat_T) 

FeaturePlot(sc21, features = feat_B, ncol=2)
VlnPlot(sc21, features = feat_B) 

FeaturePlot(sc21, features = feat_neut, ncol=2)
VlnPlot(sc21, features = feat_neut) 

FeaturePlot(sc21, features = feat_mast, ncol=2)
VlnPlot(sc21, features = feat_mast)

FeaturePlot(sc21, features = feat_DC, ncol=2)
VlnPlot(sc21, features = feat_DC) 

# We automate this process with the function AddModuleScore, which calculate an expression score for each cell:
sc21 <- AddModuleScore(sc21,features = list(feat_mac), name = "feat_mac") 
sc21 <- AddModuleScore(sc21,features = list(feat_T), name = "feat_T") 
sc21 <- AddModuleScore(sc21,features = list(feat_B), name = "feat_B") 
sc21 <- AddModuleScore(sc21,features = list(feat_neut), name = "feat_neut")
sc21 <- AddModuleScore(sc21,features = list(feat_mast), name = "feat_mast") 
sc21 <- AddModuleScore(sc21,features = list(feat_DC), name = "feat_DC")
sc21 <- AddModuleScore(sc21,features = list(eq.g2m.genes), name = "feat_eq.g2m.genes")
sc21 <- AddModuleScore(sc21,features = list(eq.s.genes), name = "feat_eq.s.genes")

# Expression patterns
sc21 <- SetIdent(sc21, value = sc21$RNA_snn_res.1)

(FeaturePlot(sc21, features = "feat_mac1")+ labs(title=NULL, x=NULL, y=NULL)) +
  (VlnPlot(sc21, features = "feat_mac1")+ labs(title=NULL, x=NULL)) &
  plot_annotation(title="Monocytes/Macrophages feature score")
# Gene expression pattern of clusters 7, 8, 9 and 15 is consistent with Mo/Ma

(FeaturePlot(sc21, features = "feat_T1")+ labs(title=NULL, x=NULL, y=NULL)) +
  (VlnPlot(sc21, features = "feat_T1")+ labs(title=NULL, x=NULL)) &
  plot_annotation(title="T cells feature score")
# Gene expression pattern of clusters 0, 1, 2, 3, 4, 5, 12 and 14 is consistent with T cells

(FeaturePlot(sc21, features = "feat_B1")+ labs(title=NULL, x=NULL, y=NULL)) +
  (VlnPlot(sc21, features = "feat_B1")+ labs(title=NULL, x=NULL)) &
  plot_annotation(title="B/Plasma cells feature score")
# Gene expression pattern of cluster 11 is consistent with B/Plasma cells

(FeaturePlot(sc21, features = "feat_neut1")+ labs(title=NULL, x=NULL, y=NULL)) +
  (VlnPlot(sc21, features = "feat_neut1")+ labs(title=NULL, x=NULL)) &
  plot_annotation(title="Neutrophils feature score")
# Gene expression pattern of cluster 6 is consistent with neutrophils

(FeaturePlot(sc21, features = "feat_mast1")+ labs(title=NULL, x=NULL, y=NULL)) +
  (VlnPlot(sc21, features = "feat_mast1")+ labs(title=NULL, x=NULL)) &
  plot_annotation(title="Mast cells feature score")
# Gene expression pattern of cluster 10 is consistent with mast cells

(FeaturePlot(sc21, features = "feat_DC1")+ labs(title=NULL, x=NULL, y=NULL)) +
  (VlnPlot(sc21, features = "feat_DC1")+ labs(title=NULL, x=NULL)) &
  plot_annotation(title="Dendritic cells feature score")
# Gene expression pattern of cluster 13 is consistent with dendritic cells

(FeaturePlot(sc21, features = "feat_eq.g2m.genes1")+ labs(title=NULL, x=NULL, y=NULL)) +
  (VlnPlot(sc21, features = "feat_eq.g2m.genes1")+ labs(title=NULL, x=NULL)) &
  plot_annotation(title="G2M phase feature score")
# Gene expression pattern of clusters 14 and 15 is consistent with proliferating cells


# MARKER GENES FOR CELL CLUSTERS

# Find the markers for each cluster at the selected resolution
sc21 <- SetIdent(sc21, value = sc21$RNA_snn_res.1)
markers_all <- FindAllMarkers(sc21, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all <- subset(markers_all, markers_all$p_val_adj < 0.05) #filtering the non significant marker genes
# Supplementary Table 2:
write.csv(markers_all, 'sc21_markers_all.csv') 

# Merge the cell clusters into the major cell types previously identified
sc21 <- RenameIdents(object=sc21,'0'="T cells",'1'= "T cells",'2'="T cells",'3'="T cells",
                                '4'="T cells",'5'="T cells",'12'="T cells",'14'="T cells",
                                '7'="Mo/Ma", '8'="Mo/Ma", '9'="Mo/Ma",'15'="Mo/Ma",
                                '6'="Neutrophils", '10'="Mast cells",'11'="B/Plasma cells",'13'="Dendritic cells" )

# Store the current identities in a new column of meta.data called CellType
sc21$CellType <- Idents(object = sc21)

# Find the markers for the major cell types 
sc21_major <- sc21
sc21_major <- SetIdent(sc21_major, value = sc21_major$CellType)
markers_all_major <- FindAllMarkers(sc21_major, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all_major <- subset(markers_all_major, markers_all_major$p_val_adj < 0.05) #filtering the non significant genes
# Supplementary Table 2:
write.csv(markers_all_major, 'sc21_markers_all_major.csv')


# QUALITY CONTROL
# Exploration of ribosomal protein gene and mitochondrial gene expression in the major cell groups
sc21_major<- PercentageFeatureSet(sc21_major, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")
VlnPlot(sc21, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeaturePlot(sc21_major, features = "percent.mito") + VlnPlot(sc21_major, features = "percent.mito")
# Supplementary Figure 2:
tiff(file='sc21_RP genes.tiff', width = dpi*12, height = dpi*5, units = "px",res = dpi)
(FeaturePlot(sc21_major, features = "percent.ribo") + VlnPlot(sc21_major, features = "percent.ribo") +
    plot_layout(ncol = 2)) & 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), legend.title = element_blank()) &
  labs(title="") &
  plot_annotation(tag_levels = 'A')
dev.off()

# FIGURE 3

# Figure 3A: UMAP with the major cell populations annotated
# create custom function to label clusters
GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
    y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
  } else {
    x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
    y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
  }
  return(list('x' = x, 'y' = y))
}

custom.LabelClusters <- function(
  plot, # Use DimPlot to generate base ggplot to apply function
  id,   # The seurat cluster identifier
  clusters = NULL,
  labels = NULL,
  split.by = NULL,
  repel = F,
  colors = colors,
  circle.size = circle.size,
  text.size = text.size,
  ...
) {
  xynames <- unlist(x = GetXYAesthetics(plot = plot), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = 'rbind',
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      return(data.medians)
    }
  )
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  geom.use <- ifelse(test = repel, yes = geom_point, no = geom_label)
  plot <- plot + geom.use(
    data = labels.loc, size = circle.size, shape=21, stroke = 0.66, col = "gray17",
    mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
    ...
  ) + geom_text( 
    size = text.size,
    data = labels.loc, col = "gray17",
    mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id),
    ...
  )
  return(plot)
}

# plotting the UMAP
colors3<-c("#0080FE", "#78BFEA","#00BED3","#5BC1BF","#678297","#489ECE","#9675B4",
           "#692E14","#F57D20","#FFCD03","#DD1A21","#F05729","#84C8E2","#9ACA3A","#5BC1BF","#F57D20")
sc21_major<-SetIdent(sc21_major,value = sc21_major$RNA_snn_res.1)
p_sc21 <-DimPlot(object = sc21_major,pt.size = 0.01,label = F,group.by = "RNA_snn_res.1") +
  scale_color_manual(values = colors3) + scale_fill_manual(values = colors3) + 
  xlab(expression('UMAP'[1])) + ylab(expression('UMAP'[2])) +
  theme(legend.position = "none") + labs(title="")
p_sc21 <- custom.LabelClusters(plot = p_sc21, 
                           id = "RNA_snn_res.1", 
                           clusters = levels(sc21_major@meta.data[,"RNA_snn_res.1"]),
                           circle.size = 5, 
                           text.color = "black",
                           text.size = 3,
                           shape = 21,
                           fill = colors3,
                           repel = T)
dpi=300
tiff(file='sc21_UMAP_major.tiff', width = dpi*5, height = dpi*4, units = "px",res = dpi)
print(p_sc21)
dev.off()

# Figure 3B: Dot plot of the 5 most upregulated genes in each major cell population.
major_top5 <- markers_all_major %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_all_major), value = TRUE)), n = 5)
major_top5$gene <- recode(major_top5$gene,"LOC102148710"="*WFDC2","LOC100060608"="*IGLC7","LOC102147726"="*IGLC1")
markers <- c("*WFDC2", "*IGLC7","*IGLC1", "RET")
markers.loc <- c("LOC102148710","LOC100060608","LOC102147726")
sc21_fig <- sc21_major
sc21_fig@meta.data[,markers]<-GetAssayData(sc21_fig)[markers.loc,] %>% as.matrix %>% t
dpi=300
tiff(file='sc21_major_dotplot.tiff', width = dpi*4, height = dpi*8, units = "px",res = dpi)
DotPlot(sc21_fig, features = unique(major_top5$gene)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(size = 10), legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
  coord_flip()
dev.off()

# Figure 3C: Gene expression patterns used for cell type assignment
dpi=300
tiff(file='sc21_featureplots.tiff', width = dpi*12, height = dpi*7, units = "px",res = dpi)
(((FeaturePlot(sc21, features = "feat_T1")+ labs(title="T cell feature expression score \n (CD2, CD3D, CD3E, CD3G)"))+
    (FeaturePlot(sc21, features = "feat_mac1")+ labs(title="Mo/Ma feature expression score \n (CD68, CD163)"))+
    (FeaturePlot(sc21, features = "feat_neut1")+ labs(title="Neutrophil feature expression score \n (CSF3R, LILRA5, RGS2, TG)"))+
    (FeaturePlot(sc21, features = "feat_mast1")+ labs(title="Mast cell feature expression score \n (GCSAML, HPGDS, LTC4S,MS4A2)"))+
    (FeaturePlot(sc21, features = "feat_B1")+ labs(title="B/Plasma cell feature expression score \n (CD79A, CD79B, MS4A1)"))+
    (FeaturePlot(sc21, features = "feat_DC1")+ labs(title="Dendritic cell feature expression score \n (CCR7, CD83)")))+
    plot_layout(ncol = 3)) &
  theme(plot.title=element_text(size=10), axis.title = element_text(size=10))
dev.off()

# Figure 3D: bar chart showing the distributions of the major bronchoalveolar cell populations obtained with cytology on BALF, 
# with cytology on the cell suspension (post cryopreservation) and with scRNA-seq on the cell suspension. 
# T cells and B/Plasma cells are counted as lymphocytes, while Mo/Ma and DCs are counted as macrophages.
sc21_major.DCC = sc21_major
sc21_major.DCC <- RenameIdents(object = sc21_major.DCC, 'B/Plasma cells'="Lymphocytes",'T cells'="Lymphocytes",
                               'Dendritic cells'="Macrophages",'Mo/Ma'="Macrophages")
prop.table(table(Idents(sc21_major.DCC), sc21_major.DCC$orig.ident),2) #values are input into DCC_pilot.txt (type = sc)
dcc.Table<-read.table("DCC_pilot.txt",header=T,sep="\t")
dcc.Table<-dcc.Table[,c(2:7)]
dcc.Table<-as.data.frame(dcc.Table)
dcc.Table.melt<-melt(dcc.Table, id=c("horse","type"))
dcc.Table.melt$horse <- as.character(dcc.Table.melt$horse)
dcc.Table.melt.group<-dcc.Table.melt %>% group_by(type)

dpi=300
tiff(file='sc21_major_DCC.tiff', width = dpi*8, height = dpi*4, units = "px",res = dpi)
ggplot(dcc.Table.melt.group,aes(x=type,y=value,fill=variable,group=type)) +
  geom_col(position = "fill", width = 0.5) +
  theme_bw(base_size = 15) +
  facet_wrap(~horse,nrow=6) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(breaks=c("cytopost", "cytopre", "sc"),
                   labels=c("BALF cytology", "Cell suspension cytology", "scRNA-seq")) +
  theme(axis.title.x = element_text(size=11), legend.title = element_blank()) + coord_flip() +
  scale_fill_manual(name="Cell type", values=c('#3399FF', '#FFCC33','#9966FF','#CC0000'),
                    labels=c("Lymphocytes", "Macrophages", "Neutrophils", "Mast cells")) +
  labs(x=NULL, y="Fraction of cells")
dev.off()

# Figure 3E: bar chart showing the distribution of the major bronchoalveolar cell populations for each horse. 
sc21_major <- SetIdent(sc21_major, value = sc21_major$CellType)
pt <- table(Idents(sc21_major), sc21_major$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt <- pt %>%
  mutate(Var2 = recode(Var2, ss1 = 'Horse 1', ss2 = 'Horse 2', ss3 =  'Horse 3' ))
pt$Var2 <- factor(pt$Var2 ,levels = c("Horse 3", "Horse 2", "Horse 1"))
tiff(file='sc21_major_DCC_orig.tiff', width = dpi*8, height = dpi*4, units = "px",res = dpi)
ggplot(pt, aes(x = Var2, y = Freq, fill = factor(Var1, levels=c("Dendritic cells","B/Plasma cells","Mast cells","Neutrophils","Mo/Ma","T cells")))) +
  labs(x=NULL, y="Fraction of cells", fill=NULL)+
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = c('#FF9900','#66CC33','#CC0000','#9966FF','#FFCC33','#3399FF')) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_text(size=11)) +
  coord_flip()
dev.off()



# INDEPENDENT ANALYSIS OF SELECTED MAJOR CELL GROUPS (SUBCLUSTERING)

# We reanalyze independently the cell groups with a high cell number and/or expected to be constituted of distinct cell subtypes:
# monocytes/macrophages, T cells and B/plasma cells.
# For each cell subset, we rerun the following steps: scaling, variable features selection, PCA, normalization and clustering. 


# MONOCYTES-MACROPHAGES (Mo/Ma)

sc21 <- SetIdent(sc21, value = sc21$CellType)
sc21_moma <- subset(sc21, idents = "Mo/Ma") 

# Normalization
sc21_moma <- NormalizeData(sc21_moma, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features selection
sc21_moma <- FindVariableFeatures(sc21_moma, selection.method = "vst", nfeatures = 2000)
vf_top10_moma <- head(VariableFeatures(sc21_moma), 10)
vf_plot_moma <- VariableFeaturePlot(sc21_moma)
LabelPoints(plot = vf_plot_moma, points = vf_top10_moma, repel = TRUE, xnudge = 0, ynudge = 0)

# Scaling
sc21_moma <- ScaleData(sc21_moma)

# Dimensionality reduction
sc21_moma <- RunPCA(sc21_moma)
DimPlot(sc21_moma, reduction = "pca")
DimPlot(sc21_moma, reduction = "pca", group.by = "Phase")  
# Cells are not clustering based on cell cycle phase; thus regression is not necessary.
# Selection of number of PCs...
#...visually with an elbow plot
ElbowPlot(sc21_moma, ndims = 40) 
#...or with a quantitative technique
pct.moma <- sc21_moma[["pca"]]@stdev / sum(sc21_moma[["pca"]]@stdev)*100
cumu.moma <- cumsum(pct.moma)
co1.moma <- which(cumu.moma > 90 & pct.moma < 5)[1]
co1.moma #43
co2.moma <- sort(which((pct.moma[1:length(pct.moma) - 1] - pct.moma[2:length(pct.moma)]) > 0.1), decreasing = T)[1] + 1
co2.moma #13
pcs.moma <- min(co1.moma, co2.moma)
plot_df.moma <- data.frame(pct.moma = pct.moma, 
                              cumu.moma = cumu.moma, 
                              rank = 1:length(pct.moma))
ggplot(plot_df.moma, aes(cumu.moma, pct.moma, label = rank, color = rank > pcs.moma)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct.moma[pct.moma > 5]), color = "grey") +
  theme_bw() + ggtitle("Elbow plot (quantitative)")
# We choose 13 PCs
DimHeatmap(sc21_moma, dims = 1:13, cells = 500, balanced = TRUE)
VizDimLoadings(sc21_moma, dims = 1, ncol = 1) + theme_minimal(base_size = 8) 

# Clusters visualization with UMAP
sc21_moma <- RunUMAP(sc21_moma, dims = 1:13) 
DimPlot(sc21_moma, reduction = "umap")
DimPlot(sc21_moma, reduction = "umap", group.by = "Phase")

# clustering 
sc21_moma <- FindNeighbors(sc21_moma, dims = 1:13) 
sc21_moma <- FindClusters(sc21_moma, resolution = seq(0.1, 1.4, by=0.1))
head(sc21_moma@meta.data)
# Visualize how clusters sub-divide at increasing resolution:
clustree(sc21_moma@meta.data[,grep("RNA_snn_res", colnames(sc21_moma@meta.data))], prefix = "RNA_snn_res.")
(DimPlot(object = sc21_moma, group.by=grep("res",colnames(sc21_moma@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Mo/Ma clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc21_moma, group.by=grep("res",colnames(sc21_moma@meta.data),value = TRUE)[9:12], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Mo/Ma clustering resolutions')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
# We choose resolution = 1.2.
DimPlot(sc21_moma, reduction = "umap", group.by = "RNA_snn_res.1.2", label=T)

# Monocytes/Macrophages cell subtypes annotation
sc21_moma <- SetIdent(sc21_moma, value = sc21_moma$RNA_snn_res.1.2) 
sc21_moma <- RenameIdents(sc21_moma, "0"="Mo/Ma 0", "1"="Mo/Ma 1", "2"="Mo/Ma 2", "3"="Mo/Ma 3", 
                         "4"="Mo/Ma 4", "5"="Mo/Ma 5", "6"="Mo/Ma 6")
sc21_moma$CellSubtype <- Idents(object = sc21_moma)

# Find the markers for the Mo/Ma cell subtypes 
sc21_moma <- SetIdent(sc21_moma, value = sc21_moma$CellSubtype) 
markers_all_moma <- FindAllMarkers(sc21_moma, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all_moma <- subset(markers_all_moma, markers_all_moma$p_val_adj < 0.05) #filtering the non significant genes
# Supplementary Table 4
write.csv(markers_all_moma, 'moma_markers_all.csv')

# Exploration of ribosomal protein gene and mitochondrial gene expression
sc21_moma<- PercentageFeatureSet(sc21_moma, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")
FeaturePlot(sc21_moma, features = "percent.ribo") + VlnPlot(sc21_moma, features = "percent.ribo")
FeaturePlot(sc21_moma, features = "percent.mito") + VlnPlot(sc21_moma, features = "percent.mito")

# The low feature and RNA count in MoMa 6 suggest there are mostly dead or dying cells
VlnPlot(sc21_moma, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Create vectors of marker genes (features) for cell subtypes 
feat_AM <- c("MARCO", "APOE", "MSR1", "CD163" )
feat_HLA <- c("DQB", "DQB.1", "DQA", "DQA.1","DRA", "DRB") 
feat_cp <- c("C4BPA", "C1QB", "C1QA", "C1QC") #complement
feat_T <- c("CD2", "CD3D", "CD3E", "CD3G")

# We automate this process with the function AddModuleScore, which calculate an expression score for each cell:
sc21_moma <- AddModuleScore(sc21_moma,features = feat_AM, name = "feat_AM")
sc21_moma <- AddModuleScore(sc21_moma,features = feat_HLA, name = "feat_HLA")
sc21_moma <- AddModuleScore(sc21_moma,features = feat_cp, name = "feat_cp")
sc21_moma <- AddModuleScore(sc21_moma,features = feat_T, name = "feat_T")
sc21_moma <- AddModuleScore(sc21_moma,features = list(eq.g2m.genes), name = "feat_eq.g2m.genes")
sc21_moma <- AddModuleScore(sc21_moma,features = list(eq.s.genes), name = "feat_eq.s.genes")
# We create a filtered object without MoMa6 (dead cells)
sc21_moma_flt <- subset(sc21_moma, idents = c("Mo/Ma 0","Mo/Ma 1","Mo/Ma 2","Mo/Ma 3","Mo/Ma 4","Mo/Ma 5"))  

# Gene expression patterns
FeaturePlot(sc21_moma_flt, features = "feat_AM1")
VlnPlot(sc21_moma_flt, features = "feat_AM1")

FeaturePlot(sc21_moma_flt, features = "feat_HLA1")
VlnPlot(sc21_moma_flt, features = "feat_HLA1")

FeaturePlot(sc21_moma_flt, features = "feat_cp1")
VlnPlot(sc21_moma_flt, features = "feat_cp1")

FeaturePlot(sc21_moma_flt, features = "feat_T1")
VlnPlot(sc21_moma_flt, features = "feat_T1")
 
FeaturePlot(sc21_moma_flt, features = "LOC100051526") # annotated as CD16
VlnPlot(sc21_moma_flt, features = "LOC100051526")
  
FeaturePlot(sc21_moma_flt, features = "feat_eq.g2m.genes1")
VlnPlot(sc21_moma_flt, features = "feat_eq.g2m.genes1")

#FIGURE 4

# Figure 4A: UMAP visualization of the 6 Mo/Ma clusters identified.
sc21_moma<-SetIdent(sc21_moma,value=sc21_moma$RNA_snn_res.1.2)
moma_colors<-c("#ED7014","#FAC80A","#80400B","#EC9706","#FCAE1E","#D23B05","#DEB887")
p_moma <-DimPlot(object = sc21_moma,pt.size = 0.05,label = TRUE,group.by = "RNA_snn_res.1.2") + 
  scale_color_manual(values = moma_colors) + scale_fill_manual(values = moma_colors) + 
  xlab(expression('UMAP'[1])) + ylab(expression('UMAP'[2])) + 
  theme(legend.position = "none") + labs(title="")
p_moma <- custom.LabelClusters(plot = p_moma, 
                           id = "RNA_snn_res.1.2", 
                           clusters = levels(sc21_moma@meta.data[,"RNA_snn_res.1.2"]),
                           circle.size = 5, 
                           text.color = "black",
                           text.size = 3, 
                           shape = 21,
                           fill = moma_colors,
                           repel = T)
dpi=300
tiff(file='moma_umap.tiff', width = dpi*5, height = dpi*4, units = "px",res = dpi)
print(p_moma)
dev.off()

# Figure 4B: violin plot of the RNA count per Mo/Ma cluster. 
dpi=300
tiff(file='moma_RNA.tiff', width = dpi*10, height = dpi*4, units = "px",res = dpi)
VlnPlot(sc21_moma, features = c("nCount_RNA"), ncol = 2) + labs(title="RNA count") +
  theme(axis.title=element_blank(), plot.title=element_text(size=10), legend.position = 'none') +
  scale_fill_brewer(palette="Purples")
dev.off()

# Figure 4C: gene expression patterns used for cell subtype assignment
dpi=300
tiff(file='moma_featureplots.tiff', width = dpi*8, height = dpi*8, units = "px",res = dpi)
(((FeaturePlot(sc21_moma, features = "feat_AM1")+ labs(title="Alveolar macrophage feature expression score \n (APOE, CD163, MARCO, MSR1)"))+
     (FeaturePlot(sc21_moma, features = "LOC100051526")+ labs(title="*CD16"))+
    (FeaturePlot(sc21_moma, features = "feat_eq.g2m.genes1")+ labs(title="G2M genes")))+
    (FeaturePlot(sc21_moma, features = "feat_HLA1")+ labs(title="MHC II feature expression score \n (DQB, DQB.1, DQA, DQA.1, DRA, DRB)"))+
        (FeaturePlot(sc21_moma, features = "feat_T1")+ labs(title="T cell feature expression score \n (CD2, CD3D, CD3E, CD3G)"))+
    (FeaturePlot(sc21_moma, features = "LOC100630729")+ labs(title="*Immunoglobulin Kappa Variable 2-30-like")) + 
    plot_layout(ncol = 2)) &
  theme(plot.title=element_text(size=10), axis.title = element_text(size=10))
dev.off()

# Figure 4D: Dot plot of the 5 most upregulated genes in each cluster
moma_top5 <- markers_all_moma %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_all_moma), value = TRUE)), n = 5)
moma_top5$gene <- recode(moma_top5$gene,"LOC100050560"="*IFITM1","LOC100050100"="*ORM2","LOC100069029"="*FCN1",
                         "LOC100059533"="*GSTP1","LOC100067916"="*NRMAL1","LOC100146489"="*H2AZ1")
moma_markers <- c("*IFITM1","*ORM2","*FCN1","*GSTP1","*NRMAL1","*H2AZ1")
moma_markers.loc <- c("LOC100050560","LOC100050100","LOC100069029","LOC100059533","LOC100067916","LOC100146489")
sc21_moma_fig <- sc21_moma
sc21_moma_fig@meta.data[,moma_markers]<-GetAssayData(sc21_moma_fig)[moma_markers.loc,] %>% as.matrix %>% t
dpi=300
tiff(file='moma_dotplot.tiff', width = dpi*4, height = dpi*8, units = "px",res = dpi)
DotPlot(sc21_moma_fig, features = unique(moma_top5$gene), 
                idents = c("Mo/Ma 0","Mo/Ma 1","Mo/Ma 2","Mo/Ma 3","Mo/Ma 4","Mo/Ma 5")) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(size = 10), legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
  coord_flip()
dev.off()

# Figure 4E: bar chart showing the distribution of the Mo/Ma clusters for each horse
pt_moma <- table(Idents(sc21_moma_flt), sc21_moma_flt$orig.ident) 
pt_moma <- as.data.frame(pt_moma)
pt_moma$Var1 <- as.character(pt_moma$Var1)
pt_moma <- pt_moma %>%
  mutate(Var2 = recode(Var2, ss1 = 'Horse 1', ss2 = 'Horse 2', ss3 =  'Horse 3' ))
pt_moma$Var2 <- factor(pt_moma$Var2 ,levels = c("Horse 3", "Horse 2", "Horse 1"))
dpi=300
tiff(file='moma_orig.tiff', width = dpi*8, height = dpi*4, units = "px",res = dpi)
ggplot(pt_moma, aes(x = Var2, y = Freq, fill = factor(Var1, levels=c("Mo/Ma 5","Mo/Ma 4","Mo/Ma 3","Mo/Ma 2","Mo/Ma 1","Mo/Ma 0")))) +
  geom_col(position = "fill", width = 0.5) +
  labs(x=NULL, y="Fraction of cells", fill=NULL) +
  scale_fill_manual(values = brewer.pal(6, "BrBG")) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_text(size=11)) +
  coord_flip()
dev.off()

# Figure 4F: Distribution plot showing T cell and intermediate monocyte (Mo/Ma 2) gene scores in the T cell, Mo/Ma2 and Mo/Ma4 clusters. 
# T cell genes score
sc21_moma4_moma2<-subset(x=sc21_moma_flt,idents=c("Mo/Ma 4", "Mo/Ma 2"))
sc21_moma4_2_T_merge<-merge(x=sc21_moma4_moma2,y=sc21_T_flt)
my_comparisons <- list( c("T cells", "Mo/Ma 4"), c("Mo/Ma 4", "Mo/Ma 2"), c("T cells", "Mo/Ma 2") )
md <-sc21_moma4_2_T_merge@meta.data %>% as.data.table()
rmcols<-colnames(md)[52:99]
for(i in rmcols) {
  sc21_moma4_2_T_merge[[i]] <- NULL 
}
sc21_moma4_2_T_merge<-AddModuleScore(sc21_moma4_2_T_merge,features =t_genes , name = "Tcellgenes")
sc21_moma4_2_T_merge<-AddModuleScore(sc21_moma4_2_T_merge,features =moma2_genes , name = "moma2cellgenes")
sc21_moma4_2_T_merge <- RenameIdents(sc21_moma4_2_T_merge, "T0"= "T cells", "T1"="T cells", "T2"="T cells", "T3" ="T cells", "T4"="T cells","T5" ="T cells","T6" ="T cells","T8" ="T cells")
sc21_moma4_2_T_merge$CellTypeforScore <- Idents(object = sc21_moma4_2_T_merge)
sc21_moma4_2_T_merge$CellTypeforScore <- factor(
  sc21_moma4_2_T_merge$CellTypeforScore,
  levels = c("T cells", "Mo/Ma 4", "Mo/Ma 2"))
Txd<-sc21_moma4_2_T_merge@meta.data[,c(52:100)] %>% group_by(sc21_moma4_2_T_merge$CellTypeforScore) %>% summarize_all(list(mean=mean))
Txdf<-t(as.data.frame(Txd))
rownames(Txdf)<-NULL
colnames(Txdf) <- Txdf[1,]
Txdf <- Txdf[-1, ] 
Txdf<-as.data.frame(Txdf)
Txdf$`T cells`<-as.numeric(Txdf$`T cells`)
Txdf$`Mo/Ma 4`<-as.numeric(Txdf$`Mo/Ma 4`)
Txdf$`Mo/Ma 2`<-as.numeric(Txdf$`Mo/Ma 2`)
Txdf.melt<-reshape2::melt(Txdf)
colnames(Txdf.melt)<-c("CellType","score")

plot_Tscore <- ggplot(Txdf.melt, aes(x = CellType, y = score,color=CellType)) + 
  geom_jitter(width=0.1,shape=19) + scale_color_manual(values=c('#3399FF','#CC7722','#FABD02')) +
  stat_summary(fun = "mean", geom = "crossbar",width = .2) + 
  stat_compare_means(comparisons = my_comparisons) + 
  xlab("") + ylab("T cell genes score") + theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=11, colour="black"), 
        axis.text.y = element_text(size=11, colour="black"),
        axis.title.y = element_text(size=11, colour="black")) 

# moma_2 (intermediate monocytes) genes score
Mxd<-sc21_moma4_2_T_merge@meta.data[,c(101:149)] %>% group_by(sc21_moma4_2_T_merge$CellTypeforScore) %>% 
  summarize_all(list(mean=mean))
Mxdf<-t(as.data.frame(Mxd))
rownames(Mxdf)<-NULL
colnames(Mxdf) <- Mxdf[1,]
Mxdf <- Mxdf[-1, ] 
Mxdf<-as.data.frame(Mxdf)
Mxdf$`T cells`<-as.numeric(Mxdf$`T cells`)
Mxdf$`Mo/Ma 4`<-as.numeric(Mxdf$`Mo/Ma 4`)
Mxdf$`Mo/Ma 2`<-as.numeric(Mxdf$`Mo/Ma 2`)
Mxdf.melt<-reshape2::melt(Mxdf)
colnames(Mxdf.melt)<-c("CellType","score")
plot_monoscore <- ggplot(Mxdf.melt, aes(x = CellType, y = score,color=CellType)) + 
  geom_jitter(width=0.1,shape=19) +
  scale_color_manual(values=c('#3399FF','#CC7722','#FABD02')) +
  stat_summary(fun = "mean", geom = "crossbar",width = .2) +
  stat_compare_means(comparisons = my_comparisons) + 
  xlab("") + ylab("Mo/Ma 2 genes score") + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none",
        axis.text.x = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.title.y = element_text(size=11, colour="black")) 
dpi=300
tiff(file='moma4_geneScore.tiff', width = dpi*8, height = dpi*4, units = "px",res = dpi)
plot_Tscore + plot_monoscore
dev.off()

# Figure 4G: Heatmap representing single cell gene expression of the 10  most upregulated genes in T cells, intermediate monocytes (Mo/Ma2) and presumptive monocyte-lymphocyte complexes (Mo/Ma 4).
sc21_moma2_4_flt_downSample<-subset(x=sc21_moma_flt,idents=c("Mo/Ma 2","Mo/Ma 4"),downsample=50)
sc21_T_flt_2_downSample<-subset(x=sc21_T_flt,downsample=10)
sc21_moma2_4_T_flt_downSample<-merge(x=sc21_moma2_4_flt_downSample,y=sc21_T_flt_2_downSample)
sc21_moma2_4_T_flt_downSample <- RenameIdents(sc21_moma2_4_T_flt_downSample, "T0"= "T cells", "T1"="T cells", "T2"="T cells", "T3" ="T cells", "T4"="T cells","T5" ="T cells","T6" ="T cells","T8" ="T cells")
sc21_moma2_4_T_flt_downSample$CellTypeHeatMap <- Idents(object = sc21_moma2_4_T_flt_downSample)
sc21_moma2_4_T_flt_downSample$CellTypeHeatMap <- factor(
  sc21_moma2_4_T_flt_downSample$CellTypeHeatMap,
  levels = c("T cells", "Mo/Ma 4", "Mo/Ma 2")
)
t_genes<-c("CCL5","GZMA","CD3E","CTSW","CD3G","LOC100147522","PTPRCAP","CD2","ETS1","HOPX","LOC102149605","LOC111771973",
           "GZMK","CD3D","LCK","LOC100147160","CGA","FYB1","CD81","KLRK1","PRF1","LTB","KLRD1","LOC100146699","LOC111775205",
           "RRP1B","LOC102147555","SH2D1A","ITGB1","EEF1B2","LAG3","S100A4","CD8A","RPLP1","CD7","CRIP1","CD5","MYCBP2",
           "SKAP1","CTSC","LOC100630090","RORA","PSIP1","CD69","GPR25","DQB","NCR3","BCL11B","HCST")
moma2_genes<-c("SPP1","TIMP1","LYZ","IL1B","C1H15orf48","CCL15","MMP9","BTG1","CD1A4","IL32","LOC100053579","CD81","SLC2A3",
               "CD44","CD99","ELK3","KLF2","S100A6","S100A4","METRNL","C10H19orf33","LOC111775970","LOC102148267","ATP1B1",
               "DRA","ATP1B3","NFE2L2","HIF1A","SIRPA","LOC100053223","MATK","RBPJ","PLD3","NINJ1","MAFB","LIMD1","SLC25A6",
               "PTAFR","PXDC1","S100A11","SULF2","FRMD8","LOC111775967","MS4A7","RPS12","DRB","CD74","RPLP0","SGK1")
moma_4_genes<-c("CST3","DQB.1","DQB","DQA","DRA","DRB","CD74","DQA.1","LOC100630244","CPVL","CCL5","DRB.2","GZMA","PTPRCAP",
                "CD3E","CTSW","CD3G","LTB","LOC100147522","CKB","DQA.2","LOC100630729","DRB.1","CD81","STK17B","PKIB","CD2",
                "TCOF1","FYB1","GPR33","LOC111771973","BASP1","C10H19orf33","RGS1","LOC102149605","RGS10","FLT3","NR4A2",
                "ETS1","GPR183","LOC100630497","LOC102148743","HOPX","KIAA1551","PARM1","CD3D","GPIHBP1","ID2","LOC100071401")
moma2_4_t_genes<-c(t_genes,moma2_genes,moma_4_genes)
moma2_4_t_genes<-unique(moma2_4_t_genes)

dpi=300
tiff(file='moma4_heatmap.tiff', width = dpi*6, height = dpi*3, units = "px",res = dpi)
dittoHeatmap(sc21_moma2_4_T_flt_downSample, genes = moma2_4_t_genes, metas = NULL, show_rownames=F,
             annot.by = c("CellTypeHeatMap"),
             order.by = c("CellTypeHeatMap"), legend.title="XXX",
             annotation_colors = list(CellTypeHeatMap = c('T cells' = '#3399FF',"Mo/Ma 4" = "yellow", "Mo/Ma 2" = "orange")),
             heatmap.colors = colorRampPalette(c("darkblue","darkblue","darkblue","white","firebrick4","firebrick4","firebrick4"))(200),
             cluster_cols = F,
             cluster_rows = T,
             fontsize = 7)
dev.off()


# T CELLS

sc21 <- SetIdent(sc21, value = sc21$CellType)
sc21_T <- subset(sc21, idents = "T cells") 

# Normalization
sc21_T <- NormalizeData(sc21_T, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features selection
sc21_T <- FindVariableFeatures(sc21_T, selection.method = "vst", nfeatures = 2000)
vf_top10_T <- head(VariableFeatures(sc21_T), 10)
vf_plot_T <- VariableFeaturePlot(sc21_T)
LabelPoints(plot = vf_plot_T, points = vf_top10_T, repel = TRUE, xnudge=0, ynudge=0)

# Scaling
sc21_T <- ScaleData(sc21_T)

# Dimensionality reduction
sc21_T <- RunPCA(sc21_T)
DimPlot(sc21_T, reduction = "pca")
DimPlot(sc21_T, reduction = "pca", group.by = "Phase") 
# Cells are not clustering based on cell cycle phase; thus regression is not necessary.
# Choice of adequate number of PCs...
#... using a visual approach
ElbowPlot(sc21_T, ndims = 40)  +
  geom_vline(xintercept = 11, linetype="dashed", color = "purple") +
  geom_vline(xintercept = 13, linetype="dashed", color = "blue") +
  geom_vline(xintercept = 15, linetype="dashed", color = "green")  
#...or a quantitative technique
pct.T <- sc21_T[["pca"]]@stdev / sum(sc21_T[["pca"]]@stdev)*100
cumu.T <- cumsum(pct.T)
co1.T <- which(cumu.T > 90 & pct.T < 5)[1]
co1.T #44
co2.T <- sort(which((pct.T[1:length(pct.T) - 1] - pct.T[2:length(pct.T)]) > 0.1), decreasing = T)[1] + 1
co2.T #11
pcs.T <- min(co1.T, co2.T)
plot_df.T <- data.frame(pct.T = pct.T, 
                          cumu.T = cumu.T, 
                          rank = 1:length(pct.T))
ggplot(plot_df.T, aes(cumu.T, pct.T, label = rank, color = rank > pcs.T)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct.T[pct.T > 5]), color = "grey") +
  theme_bw() + ggtitle("T - Elbow plot (quantitative)")

# We choose 11 PCs
DimHeatmap(sc21_T, dims = 1:11, cells = 500, balanced = TRUE)
VizDimLoadings(sc21_T, dims = 1, ncol = 1) + theme_minimal(base_size = 8) 

# Clusters visualization with UMAP
sc21_T <- RunUMAP(sc21_T, dims = 1:11) 
DimPlot(sc21_T, reduction = "umap")
DimPlot(sc21_T, reduction = "umap", group.by = "Phase")

# Clustering
sc21_T <- FindNeighbors(sc21_T, dims = 1:11) 
sc21_T <- FindClusters(sc21_T, resolution = seq(0.1, 1.4, by=0.1))
head(sc21_T@meta.data)
# Visualize how clusters sub-divide at increasing resolution
clustree(sc21_T@meta.data[,grep("RNA_snn_res", colnames(sc21_T@meta.data))], prefix = "RNA_snn_res.")
# Visualize the UMAP for the different resolutions
(DimPlot(object = sc21_T, group.by=grep("res",colnames(sc21_T@meta.data),value = TRUE)[1:4], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions for T cells')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc21_T, group.by=grep("res",colnames(sc21_T@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions for T cells')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc21_T, group.by=grep("res",colnames(sc21_T@meta.data),value = TRUE)[9:12], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions for T cells')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
# We choose resolution = 0.5
DimPlot(sc21_T, reduction = "umap", group.by = "RNA_snn_res.0.5", label=T) 

# T cell subtypes annotation
sc21_T <- SetIdent(sc21_T, value = sc21_T$RNA_snn_res.0.5)
sc21_T <- RenameIdents(sc21_T, "0"="T0", "1"="T1", "2"="T2", "3"="T3", 
                         "4"="T4", "5"="T5", "6"="T6", "7"="T7", "8"="T8")
sc21_T$CellSubtype <- Idents(object = sc21_T)

# Find the markers for the T cell subtypes 
sc21_T <- SetIdent(sc21_T, value = sc21_T$CellSubtype) 
markers_all_T <- FindAllMarkers(sc21_T, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all_T <- subset(markers_all_T, markers_all_T$p_val_adj < 0.05) #filtering the non significant genes
# Supplementary Table 5
write.csv(markers_all_T, 'T_markers_all.csv')

# Exploration of ribosomal protein gene and mitochondrial gene expression
sc21_T<- PercentageFeatureSet(sc21_T, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")
FeaturePlot(sc21_T, features = "percent.ribo") + VlnPlot(sc21_T, features = "percent.ribo")
FeaturePlot(sc21_T, features = "percent.mito") + VlnPlot(sc21_T, features = "percent.mito")

# The low feature and RNA count in MoMa 6 suggest there are mostly dead or dying cells
VlnPlot(sc21_T, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# We create a filtered object without T 7 (dead cells)
sc21_T_flt <- subset(sc21_T, idents = c("T0","T1","T2","T3","T4","T5","T6", "T8")) 

# Gene expression patterns
sc21_T <- AddModuleScore(sc21_T,features = list(eq.g2m.genes), name = "feat_eq.g2m.genes")
sc21_T <- AddModuleScore(sc21_T,features = list(eq.s.genes), name = "feat_eq.s.genes")

FeaturePlot(sc21_T_flt, features = "CD4")
VlnPlot(sc21_T_flt, features = "CD4")

FeaturePlot(sc21_T_flt, features = c("CD8A","CD8B")) 
VlnPlot(sc21_T_flt, features = c("CD8A","CD8B")) 

FeaturePlot(sc21_T_flt, features = c("FOXP3","CTLA4")) 
VlnPlot(sc21_T_flt, features = c("FOXP3","CTLA4")) 

FeaturePlot(sc21_T_flt, features = c("GZMA","GZMK", "LOC100147522","GZMM")) #LOC100147522 annotated as GZMH
VlnPlot(sc21_T_flt, features = c("GZMA","GZMK", "LOC100147522")) 
  
FeaturePlot(sc21_T_flt, features = c("PRF1","CTSW"))
VlnPlot(sc21_T_flt, features = c("PRF1","CTSW"))

FeaturePlot(sc21_T_flt, features = "LOC100066851") #LOC100066851 annotated as SCART1
VlnPlot(sc21_T_flt, features = "LOC100066851")

FeaturePlot(sc21_T_flt, features = "ITGAE")
VlnPlot(sc21_T_flt, features = "ITGAE")

FeaturePlot(sc21_T_flt, features = "CD40LG")
VlnPlot(sc21_T_flt, features = "CD40LG")

VlnPlot(sc21_T, features = c("IFNG","IL4", "IL5","IL13","IL17A","IL17B","IL17C","IL17F","IL22","IL21"))

VlnPlot(sc21_T_flt, "feat_eq.g2m.genes1") # not correlated to RP gene expression
FeaturePlot(sc21_T_flt, "feat_eq.g2m.genes1") 


# FIGURE 5
# Figure 5A: UMAP visualization of the 9 T cell clusters identified
sc21_T<-SetIdent(sc21_T,value=sc21_T$RNA_snn_res.0.5)
t_colors<-c(colors3[1:6],"#016064",colors3[15],"#0492C2")
p_T<-DimPlot(object = sc21_T,pt.size = 0.05,label = TRUE,group.by = "RNA_snn_res.0.5") + 
  scale_color_manual(values = t_colors) + scale_fill_manual(values = t_colors) + 
  xlab(expression('UMAP'[1])) + ylab(expression('UMAP'[2])) + 
  theme(legend.position = "none") + labs(title="")

p_T <- custom.LabelClusters(plot = p_T, 
                           id = "RNA_snn_res.0.5", 
                           clusters = levels(sc21_T@meta.data[,"RNA_snn_res.0.5"]),
                           circle.size = 5, 
                           text.color = "black",
                           text.size = 3, 
                           shape = 21,
                           fill = t_colors,
                           repel = T)
dpi=300
tiff(file='T_umap.tiff', width = dpi*5, height = dpi*4, units = "px",res = dpi)
print(p_T)
dev.off()

# Figure 5B: Violin plot of the RNA count per T cell cluster. 
dpi=300
tiff(file='T_RNA.tiff', width = dpi*10, height = dpi*4, units = "px",res = dpi)
VlnPlot(sc21_T, features = c("nCount_RNA"), ncol = 2) + labs(title="RNA count") +
  theme(axis.title=element_blank(), plot.title=element_text(size=10), legend.position = 'none') +
  scale_fill_brewer(palette="Purples")
dev.off()

# Figure 5C: gene expression patterns used for cell subtype assignment
dpi=300
tiff(file='T_featureplots.tiff', width = dpi*8, height = dpi*8, units = "px",res = dpi)
(((FeaturePlot(sc21_T, features = "feat_eq.g2m.genes1")+ labs(title="G2M genes"))+
    (FeaturePlot(sc21_T, features = "TCF7")+ labs(title="TCF7"))+
    (FeaturePlot(sc21_T, features = "GZMA")+ labs(title="GZMA")))+
    (FeaturePlot(sc21_T, features = "LOC100062846")+ labs(title="*KLRC1"))+
    (FeaturePlot(sc21_T, features = "LOC100066851")+ labs(title="*SCART1"))+
    (FeaturePlot(sc21_T, features = "FOXP3")+ labs(title="FOXP3")) + 
    plot_layout(ncol = 2)) &
  theme(plot.title=element_text(size=10), axis.title = element_text(size=10))
dev.off()

# Figure 5D: dot plot of the 5 most upregulated genes in each cluster 
T_top5 <- markers_all_T %>% group_by(cluster) %>%
  slice_max(get(grep("^avg_log", colnames(markers_all_T), value = TRUE)), n = 5)
T_top5$gene <- recode(T_top5$gene,"LOC100066851"="*SCART1","LOC100062823"="*KLRC1","LOC101910264"="*KLRD1",
                      "LOC100062846"="*KLRC1","LOC100051986"="*GZMH","LOC100147522"="*GZMH", "SEPT11"="SEPTIN11")
T_top5.flt <- T_top5$gene[! T_top5$gene %in% c("LOC111771923","RPL5","RPS12")]# remove rRNA and ribosomal protein genes 
T_markers <- c("*SCART1","*KLRC1","*KLRD1","*KLRC1","*GZMH","*GZMH","SEPTIN11")
T_markers.loc <- c("LOC100066851","LOC100062823","LOC101910264","LOC100062846","LOC100051986","LOC100147522","SEPT11")
sc21_T_fig <- sc21_T
sc21_T_fig@meta.data[,T_markers]<-GetAssayData(sc21_T_fig)[T_markers.loc,] %>% as.matrix %>% t
dpi=300
tiff(file='T_dotplot.tiff', width = dpi*4, height = dpi*8, units = "px",res = dpi)
DotPlot(sc21_T_fig, features = unique(T_top5.flt), 
        idents = c("T0","T1","T2","T3","T4","T5","T6","T8")) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(size = 10), legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
  coord_flip()
dev.off()

# Figure 5E:  Distribution of the T cell clusters for each horse
pt_T <- table(Idents(sc21_T_flt), sc21_T_flt$orig.ident) 
pt_T <- as.data.frame(pt_T)
pt_T$Var1 <- as.character(pt_T$Var1)
pt_T <- pt_T %>%
  mutate(Var2 = recode(Var2, ss1 = 'Horse 1', ss2 = 'Horse 2', ss3 =  'Horse 3' ))
pt_T$Var2 <- factor(pt_T$Var2 ,levels = c("Horse 3", "Horse 2", "Horse 1"))
dpi=300
tiff(file='T_orig.tiff', width = dpi*8, height = dpi*4, units = "px",res = dpi)
ggplot(pt_T, aes(x = Var2, y = Freq, fill = factor(Var1, levels=c("T8","T6","T5","T4","T3","T2","T1","T0")))) +
  geom_col(position = "fill", width = 0.5) +
  labs(x=NULL, y="Fraction of cells", fill=NULL) +
  scale_fill_manual(values = brewer.pal(8, "RdBu")) +
  scale_y_continuous(labels = scales::percent) +
   guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_text(size=11)) +
  coord_flip()
dev.off()

# Figure 5F: Violin plots showing expression of CD4, CD8A and CD8B among T cell clusters 
dpi=300
tiff(file='T_violinplots.tiff', width = dpi*12, height = dpi*3, units = "px",res = dpi)
VlnPlot(sc21_T, features = c("CD4","CD8A","CD8B")) &
  theme(axis.title=element_blank(), plot.title=element_text(size=10), legend.position = 'none') &
  scale_fill_brewer(palette="Purples", direction = "-1")
dev.off()


#B CELLS

sc21 <- SetIdent(sc21, value = sc21$CellType)
sc21_B <- subset(sc21, idents = "B/Plasma cells") 

# Normalization
sc21_B <- NormalizeData(sc21_B, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable features selection
sc21_B <- FindVariableFeatures(sc21_B, selection.method = "vst", nfeatures = 2000)
vf_top10_B <- head(VariableFeatures(sc21_B), 10)
vf_plot_B <- VariableFeaturePlot(sc21_B)
LabelPoints(plot = vf_plot_B, points = vf_top10_B, repel = TRUE, xnudge=0, ynudge=0)

# Scaling
sc21_B <- ScaleData(sc21_B)

# Dimensionality reduction
sc21_B <- RunPCA(sc21_B)
DimPlot(sc21_B, reduction = "pca")
DimPlot(sc21_B, reduction = "pca", group.by = "Phase")  
# Cells are not clustering based on cell cycle phase; thus regression is not necessary.
# Choice of adequate number of PCs...
#... using a visual approach
ElbowPlot(sc21_B, ndims = 40) 
#...or a quantitative technique
pct.B <- sc21_B[["pca"]]@stdev / sum(sc21_B[["pca"]]@stdev)*100
cumu.B <- cumsum(pct.B)
co1.B <- which(cumu.B > 90 & pct.B < 5)[1]
co1.B #44
co2.B <- sort(which((pct.B[1:length(pct.B) - 1] - pct.B[2:length(pct.B)]) > 0.1), decreasing = T)[1] + 1
co2.B #8
pcs.B <- min(co1.B, co2.B)
plot_df.B <- data.frame(pct.B = pct.B, 
                          cumu.B = cumu.B, 
                          rank = 1:length(pct.B))
ggplot(plot_df.B, aes(cumu.B, pct.B, label = rank, color = rank > pcs.B)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct.B[pct.B > 5]), color = "grey") +
  theme_bw() + ggtitle("B/Plasma cells - Elbow plot (quantitative)")
# We choose 8 PCs 
DimHeatmap(sc21_B, dims = 1:8, cells = 500, balanced = TRUE)
VizDimLoadings(sc21_B, dims = 1, ncol = 1) + theme_minimal(base_size = 8) 

# Clusters visualization with UMAP
sc21_B <- RunUMAP(sc21_B, dims = 1:8) 
DimPlot(sc21_B, reduction = "umap")
DimPlot(sc21_B, reduction = "umap", group.by = "Phase")

# Clustering
sc21_B <- FindNeighbors(sc21_B, dims = 1:8) 
sc21_B <- FindClusters(sc21_B, resolution = seq(0.1, 1.2, by=0.1))
head(sc21_B@meta.data)
# Visualize how clusters sub-divide at increasing resolution:
clustree(sc21_B@meta.data[,grep("RNA_snn_res", colnames(sc21_B@meta.data))], prefix = "RNA_snn_res.")
# Visualize the UMAP for the different resolutions
(DimPlot(object = sc21_B, group.by=grep("res",colnames(sc21_B@meta.data),value = TRUE)[5:8], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions for B/Plasma cells')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
(DimPlot(object = sc21_B, group.by=grep("res",colnames(sc21_B@meta.data),value = TRUE)[9:12], ncol=2 , pt.size=0.5, reduction = "umap", label = T) +
    plot_annotation(title = 'Clustering resolutions for B/Plasma cells')) & 
  theme(legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size=6))
# We choose resolution = 0.7
DimPlot(sc21_B, reduction = "umap", group.by = "RNA_snn_res.0.7", label=T)

# B/Plasma cell subtype annotation
sc21_B <- SetIdent(sc21_B, value = sc21_B$RNA_snn_res.0.7)
sc21_B <- RenameIdents(sc21_B, "0"="B/Plasma 0", "1"="B/Plasma 1")
sc21_B$CellSubtype <- Idents(object = sc21_B)

# Find the markers for the B/Plasma cell subtypes 
sc21_B <- SetIdent(sc21_B, value = sc21_B$CellSubtype)
markers_all_B <- FindAllMarkers(sc21_B, only.pos = F, min.pct = 0.25, thresh.use = 0.25) 
markers_all_B <- subset(markers_all_B, markers_all_B$p_val_adj < 0.05) #filtering the non significant genes
# Supplementary Table 6
write.csv(markers_all_B, 'BPlasma_markers_all.csv')

# Create vectors of marker genes (features) for cell subtypes 
feat_Bcell <- c("CD19", "CD74", "BANK1", "MS4A1", "CD79B")
feat_plasma <- c("JCHAIN") 

# Gene expression patterns
FeaturePlot(sc21_B, features=feat_Bcell) 
VlnPlot(sc21_B, features=feat_Bcell)  

FeaturePlot(sc21_B, features=feat_plasma) 
VlnPlot(sc21_B, features=feat_plasma) 

# Exploration of ribosomal protein gene and mitochondrial gene expression
VlnPlot(sc21_B, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
sc21_B<- PercentageFeatureSet(sc21_B, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", col.name = "percent.ribo")
FeaturePlot(sc21_B, features = "percent.ribo") + VlnPlot(sc21_B, features = "percent.ribo")
FeaturePlot(sc21_B, features = "percent.mito") + VlnPlot(sc21_B, features = "percent.mito")

# FIGURE 6

# Figure 6A: UMAP visualization of the 2 B/Plasma cell clusters identified
sc21_B <- SetIdent(sc21_B, value = sc21_B$CellSubtype)
b_colors<-c("#03C04A","#466D1D")
p_B <- DimPlot(object = sc21_B,pt.size = 0.05,label = TRUE,group.by = "RNA_snn_res.0.7") + 
  scale_color_manual(values = b_colors) + scale_fill_manual(values = b_colors) + 
  xlab(expression('UMAP'[1])) + ylab(expression('UMAP'[2])) + 
  theme(legend.position = "none") + labs(title="")

p_B <- custom.LabelClusters(plot = p_B, 
                           id = "RNA_snn_res.0.7", 
                           clusters = levels(sc21_B@meta.data[,"RNA_snn_res.0.7"]),
                           circle.size = 5, 
                           text.color = "black",
                           text.size = 3, 
                           shape = 21,
                           fill = b_colors,
                           repel = T)

dpi=300
tiff(file='B_umap.tiff', width = dpi*5, height = dpi*4, units = "px",res = dpi)
print(p_B)
dev.off()

# Figure 6B: bar chart showing the distribution of the B/Plasma cell clusters for each horse
sc21_B <- SetIdent(sc21_B, value = sc21_B$CellSubtype)
pt_B <- table(Idents(sc21_B), sc21_B$orig.ident) 
pt_B <- as.data.frame(pt_B)
pt_B$Var1 <- as.character(pt_B$Var1)
pt_B <- pt_B %>%
  mutate(Var2 = recode(Var2, ss1 = 'Horse 1', ss2 = 'Horse 2', ss3 =  'Horse 3' ))
pt_B$Var2 <- factor(pt_B$Var2 ,levels = c("Horse 3", "Horse 2", "Horse 1"))
dpi=300
tiff(file='BPlasma_orig.tiff', width = dpi*8, height = dpi*4, units = "px",res = dpi)
ggplot(pt_B, aes(x = Var2, y = Freq, fill = factor(Var1, levels=c("B/Plasma 1","B/Plasma 0")))) +
  geom_col(position = "fill", width = 0.5) +
  labs(x=NULL, y="Fraction of cells", fill=NULL) +
  scale_fill_manual(values = brewer.pal(3,"Greens")) +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_text(size=11)) +
  coord_flip()
dev.off()

# Figure 6C: gene expression patterns used for cell subtype assignment
dpi=300
tiff(file='B_featureplots.tiff', width = dpi*12, height = dpi*8, units = "px",res = dpi)
(((FeaturePlot(sc21_B, features = "CD79B")+ labs(title="CD79B"))+
    (FeaturePlot(sc21_B, features = "DRB")+ labs(title="DRB"))+
    (FeaturePlot(sc21_B, features = "LOC100147624")+ labs(title="*DERL3")))+
    (FeaturePlot(sc21_B, features = "JCHAIN")+ labs(title="JCHAIN"))+
    plot_layout(ncol = 2)) &
  theme(plot.title=element_text(size=15), axis.title = element_text(size=15))
dev.off()
########################################################
########################################################
# Section : clustering and DE analysis or gene markers discovery
# feature selection and dimension reduction was done already in the batch correction part
# So we start with the PCs from fastMNN in scran, with which we define distance for clustering
# 1) different clustering methods will be tested  
# 2) special design matrix will be used for DE analysis 
# http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/de.html#2_blocking_on_uninteresting_factors_of_variation
########################################################
########################################################

########################################################
# Section : Clustering section by integrating various informations: 
# gene expression, fac info, estimated timing and velocity 
########################################################
########################################################


library(scater)
library(scran)
library(scRNA.seq.funcs)
library(matrixStats)
#library(M3Drop)
library(RColorBrewer)
library(SingleCellExperiment)

load(file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_bcMNN.Rdata')) 

Seurat.clustering = TRUE
Test.scran.clustering = FALSE
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))

## double check the umap and tsnet visulization and 
set.seed(100)
sce <- runUMAP(sce, use_dimred="MNN", n_dimred = 15)
p1 = plotUMAP(sce, ncomponents = 2, colour_by="batches", size_by = "FSC_log2", point_size= 0.01) + ggtitle("Corrected-umap") 
plot(p1)

sce = runTSNE(sce, use_dimred="MNN", perplexity = 20, n_dimred = 15)
p2 = plotTSNE(sce, ncomponents = 2, colour_by="batches", size_by = "FSC_log2", point_size= 0.01) + ggtitle("Corrected-tSNE") 
plot(p2)

p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") + 
  fontsize + ggtitle("MNN corrected")
p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
  fontsize + ggtitle("MNN corrected")
multiplot(p1, p2, cols = 2)


p1 = plotTSNE(sce, colour_by="pha-4", size_by = "FSC_log2") + 
  fontsize + ggtitle("MNN corrected")
p2 = plotTSNE(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
  fontsize + ggtitle("MNN corrected")
multiplot(p1, p2, cols = 2)

##########################################
# (test clustering method) but Seurat method is used at the end 
# https://master.bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/work-1-reads.html/
# DE analysis is tested with scran but will be considered to be replaced by MAST or Seurat
# DE analysis (or marker gene finding) following the cluster analysis
##########################################
if(Seurat.clustering)
{ 
  require(Seurat)
  library(cowplot)
  
  check.individualExample.geneMarker = FALSE
  ntops = 3 # nb of top gene markers
  resolutions = c(0.4, 0.6, 0.8, seq(1.0, 4.0, by = 0.5))
  
  for(rr in resolutions){
    
    rr = 1.2
    cat("--- resolution is :", rr, "---\n")
    
    pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_batchCorrected_clustering.Seurat_geneMarkers.scran_resolution_", 
                     rr, ".pdf")
    pdf(pdfname, width=22, height = 18)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    ##########################################
    # test graph-based Louvain algorithm 
    ##########################################
    pbmc = as.Seurat(sce)
    #pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    #pbmc <- ScaleData(pbmc, features = rownames(pbmc))
    # pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
    pbmc = FindNeighbors(object = pbmc, reduction = "MNN", k.param = 20, dims = 1:20)
    pbmc = FindClusters(pbmc, resolution = rr, algorithm = 3)
    
    # note that you can set `label = TRUE` or use the LabelClusters function to help label
    # individual clusters
    #pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10, reduction = "MNN", reduction.key = "umap", n.neighbors = 15, repulsion.strength = 2)
    #DimPlot(pbmc, reduction = "umap")
    #FeaturePlot(pbmc, features = c("pha-4", "hnd-1"), reduction = "umap")
    # FeaturePlot(pbmc, features = c("pha-4", "hnd-1"), reduction = "UMAP")
    sce$cluster_seurat <- factor(pbmc@active.ident)
    sce$cluster <- factor(pbmc@active.ident)
    
    plotUMAP(sce, colour_by="cluster", size_by = "FSC_log2") + 
      fontsize + ggtitle("Seurat clustering")
    
    plotTSNE(sce, colour_by="cluster", size_by = "FSC_log2") + 
      fontsize + ggtitle("Seurat clustering")
    
    p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") + 
      fontsize + ggtitle("Seurat clustering")
    p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
      fontsize + ggtitle("Seurat clustering")
    multiplot(p1, p2, cols = 2)
    
    p1 = plotTSNE(sce, colour_by="pha-4", size_by = "FSC_log2") + 
      fontsize + ggtitle("Seurat clustering")
    p2 = plotTSNE(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
      fontsize + ggtitle("Seurat clustering")
    multiplot(p1, p2, cols = 2)
    
    my.clusters = as.numeric(as.character(sce$cluster_seurat))
    cat(table(my.clusters), "\n")
    
    ## run the find markers and then collect markers for each clusters (scran) 
    # https://bioconductor.org/packages/3.10/workflows/vignettes/simpleSingleCell/inst/doc/de.html#blocking-on-uninteresting-factors-of-variation
    #design <- model.matrix( ~ sce$batches)
    #design <- design[,-1,drop=FALSE]
    markers <- findMarkers(sce, my.clusters, block=sce$batches, direction="up")
    top.markers = c()
    
    for(n in unique(my.clusters)){
      #n = 0
      marker.set <- markers[[as.character(n)]]
      #marker.set <- markers[["1"]]
      #head(marker.set, 5)
      top.markers <- c(top.markers, rownames(marker.set)[marker.set$Top <= ntops])  
    }
    
    top.markers = unique(top.markers)
    cat(length(top.markers), " gene markers found by scran \n")
    
    plotHeatmap(sce, features=top.markers,
                columns=order(sce$cluster_seurat), 
                colour_columns_by=c("cluster",  "FSC_log2", "batches"),
                cluster_cols=FALSE, show_colnames = FALSE,
                center=TRUE, symmetric=TRUE, zlim=c(-5, 5))
    
    if(check.individualExample.geneMarker){
      for(n in 1:length(top.markers)) {
        xx = plotTSNE(sce, colour_by = top.markers[n]) 
        plot(xx)
      }
    }
    
    dev.off()
    
  }
  
}

##########################################
# here select subset of whole dataset to redo clustering
# and marker gene finding
# this is the beginning of iterative clustering, 
# the gene marker finding will be a criterion to stop iteration
##########################################
Select.Early.timePoints = FALSE
if(Select.Early.timePoints){
  
  pbmc = as.Seurat(sce)
  #Seurat::DimPlot(pbmc, dims = c(1, 2), reduction = "MNN")
  pbmc = FindNeighbors(object = pbmc, reduction = "MNN", k.param = 10, dims = 1:10)
  pbmc = FindClusters(pbmc, resolution = 2, algorithm = 3)
  sce$cluster_seurat <- factor(pbmc@active.ident)
  sce$cluster <- factor(pbmc@active.ident)
  
  xx = table(sce$cluster,sce$Batch)
  #colnames(xx) = sce$Batch
  cluster4early = rownames(xx)[which(xx[, 1]>=5|xx[,2]>=5)]
  
  mm = match(sce$cluster, factor(cluster4early))
  
  sels = which(!is.na(mm))
  
  cat("estimated cell in early stage -- ", length(cluster4early), 
      "clusters with", length(sels),  "cells\n")
  
  sce.sel = sce[, sels ]
  # here rerun the clustering for clusters or stages
  pcs <- reducedDim(sce.sel, "MNN")
  my.dist <- dist(pcs[, c(1:20)])
  my.tree <- hclust(my.dist, method="ward.D2")
  
  #hist(my.tree)
  library(dynamicTreeCut)
  my.clusters <- unname(cutreeDynamic(my.tree,
                                      distM=as.matrix(my.dist), 
                                      minClusterSize=5, verbose=0))
  
  table(my.clusters, sce.sel$Batch)
  
  sce.sel$cluster <- factor(my.clusters)
  
  set.seed(100)
  sce.sel <- runTSNE(sce.sel,  use_dimred = "MNN", perplexity = 20, n_dimred = 20)
  
  set.seed(100)
  sce.sel = runUMAP(sce.sel, use_dimred="MNN", perplexity = 20, n_dimred = 20)
  
  plotTSNE(sce.sel, colour_by="cluster", size_by = "total_features_by_counts") + 
    fontsize + ggtitle("scran -- hcluster")
  plotUMAP(sce.sel, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
    fontsize + ggtitle("scran -- hcluster")
  
  design <- model.matrix( ~ sce.sel$Batch)
  design <- design[,-1,drop=FALSE]
  
  markers <- findMarkers(sce.sel, sce.sel$cluster, design=design, direction = 'any')
  
  ntops = 5;
  top.markers = c()
  
  for(n in unique(my.clusters)){
    #n = 0
    marker.set <- markers[[as.character(n)]]
    #marker.set <- markers[["1"]]
    #head(marker.set, 5)
    top.markers <- c(top.markers, rownames(marker.set)[marker.set$Top <= ntops])  
  }
  top.markers = unique(top.markers)
  
  pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_batchCorrected_clustering_markerGenes_earlyTimepoint_tSNEexample.pdf")
  pdf(pdfname, width=12, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  plotTSNE(sce.sel, colour_by="cluster", size_by = "total_features_by_counts") + 
    fontsize + ggtitle("scran -- hcluster (tSNE)")
  plotUMAP(sce.sel, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
    fontsize + ggtitle("scran -- hcluster (UMAP)")
  
  plotHeatmap(sce.sel, features=top.markers,
              columns=order(sce.sel$cluster), 
              colour_columns_by=c("cluster"),
              cluster_cols=FALSE, show_colnames = FALSE,
              center=TRUE, symmetric=TRUE, zlim = c(-5, 5))
  
  for(n in 1:length(top.markers)) {
    xx = plotTSNE(sce.sel, colour_by = top.markers[n]) 
    plot(xx)
  }
  
  dev.off()
  
}
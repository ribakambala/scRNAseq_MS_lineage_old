########################################################
########################################################
# Section : Dealing with batch confouner
# first thing to try is the MNN method from scran pacakge after contacting with the first 
# author, Laleh Haghverdi and Aaron Lun
# !!!!!! to correct the batch, the main ideal is to 
# 1) select highly variable genes 2) using these genes to reduce the dimensions and correct the batch in the low dimension
# 3) the corrected output will be low dimensional output, which will be used to define the distance for clustering and also trajectory
# or the corrected gene expression matrix (but not used for DE analysis), from which PCs can also be found by PCA
# Taken together, HGV, dimensionality reduction, batch correction are together. 
# cluster and DE analysis will be the next two steps. 
# in addition, the DE analysis will be performed with original counts by providing cluster labels and batches(
# http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/de.html#2_blocking_on_uninteresting_factors_of_variation)
########################################################
########################################################
version.DATA = 'R6875_R7116_R7130_R7130redo_R7133_scRNA_v1'
version.analysis =  paste0(version.DATA, '_20190506')

dataDir = paste0("../data/")
resDir = paste0("../results/", version.analysis)
tabDir = paste0("../results/", version.analysis, "/tables/")
RdataDir = paste0("../results/", version.analysis, "/Rdata/")
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

correct.cellCycle = TRUE

##########################################
# Remove the cell cycle confounder 
# here we choose to use Seurat to regress out the cell cycle effect
# we need to train the cells to identify the cell cycle phase
# this could be more complicated than expected
##########################################
if(correct.cellCycle){
  source("scRNAseq_functions.R")
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))
  
  # cellCycle.correction(sce, method = "seurat")
  #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected.Rdata'))
  #sce_old = sce
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2.Rdata')) 
  #saveRDS(sce, file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected.rds'))
  
  library(scater)
  p1 = scater::plotPCA(
    sce,
    run_args = list(exprs_values = "logcounts", feature_set = c(1:500)), 
    size_by = "total_counts",
    #size_by = "total_features_by_counts",
    colour_by = "seqInfos"
  ) + ggtitle(paste0("PCA -- "))
  
  p2 = scater::plotPCA(
    sce,
    run_args = list(exprs_values = "logcounts_seurat", feature_set = c(1:500)), 
    size_by = "total_counts",
    #size_by = "total_features_by_counts",
    colour_by = "seqInfos"
  ) + ggtitle(paste0("PCA -- "))
  
  multiplot(p1, p2, cols = 2)
  
}

##########################################
# here to prepare the inputs for MNN in scran 
# Feature selection (HVGs) is performed
# here we just use the scran trandVar()
# further test to follow for other methods
# (https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby011/4898116)
# But many of them, e.g. BASics, Brennecke works better with spike-in 
##########################################
library(scRNA.seq.funcs)
library(scater)
library(SingleCellExperiment)
library(scran)
library(kBET)
set.seed(1234567)
options(stringsAsFactors = FALSE)

#load(file = file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos.Rdata'))
Use.FACS.informations = TRUE

##########################################
# Here we start to add the facs information in the metadata
# 
##########################################
if(Use.FACS.informations){
  #source("scRNAseq_functions.R")
  #sce = Integrate.FACS.Information(sce)
  #save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrectedv2_facsInfos.Rdata')) 
  load(file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos.Rdata'))
}else{
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2.Rdata')) 
}

# logcounts(sce) =  assay(sce, "logcounts_seurat_SG2MCorrected")
sce = sce[, which(sce$nb.cells == 1)]
sce$FSC_log2 = 3/2*log2(sce$FSC)
sce$BSC_log2 = 3/2*log2(sce$BSC)


plotColData(sce,
            x = "FSC_log2",
            y = "BSC_log2",
            colour_by = "seqInfos",
            shape_by = "seqInfos"
            
)

##########################################
# Feature selection (select HGVs) for batch corrections
# there are two options: batch-specific or use batch as block
# https://www.bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html
##########################################
batches = sce$seqInfos # choose the batches (either plates or request)
bc.uniq = unique(batches)
sce$batches <- batches

Use.fastMNN = TRUE
Norm.Vars.per.batch = FALSE # HVGs for each batch or not 
Rescale.Batches = TRUE # scaling data in each batch or not 
k.mnn = 20
cos.norm = TRUE
nb.pcs = 50
order2correct = c(4, 3, 1, 2)
#order2correct = c(3, 4, 1, 2)

source("scRNAseq_functions.R")
gene.chosen = find.HVGs(sce, Norm.Vars.per.batch = Norm.Vars.per.batch, method = "scran", ntop = 2000)
cat("nb of HGV : ", length(gene.chosen), "\n")

##########################################
# Batch correction using fastMNN from scran
# here we are using fastMNN
##########################################
if(Use.fastMNN){
  ## rescaling for each batch is recommended by the author
  ## We adjust the size factors with multiBatchNorm() to make them more comparable across batches. 
  ## This mitigates differences in scale and variance in the log-expression values between batches, especially between technologies.
  ## https://www.bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/multibatch.html#3_feature_selection_across_batches
  if(Rescale.Batches){
    ttxt = c("nout = multiBatchNorm(")
    for(n in 1:length(bc.uniq)){
      if(n != length(bc.uniq)) {
        ttxt = paste0(ttxt, "sce[, which(sce$batches == '", bc.uniq[n], "')], ")
      } else {
        ttxt = paste0(ttxt, "sce[, which(sce$batches == '", bc.uniq[n], "')])")
      }
    }
    eval(parse(text = ttxt)) ## rescaling done
    
    kk2check = 2
    par(mfrow=c(1,1))
    plot(sizeFactors(nout[[kk2check]]), sizeFactors(sce[, which(sce$batches == bc.uniq[kk2check])])); abline(0, 1, col='red')
    
  }
  original = list()
  fscs = c()
  #original0 = list()
  for(n in 1:length(bc.uniq)){
    #xx = nout[[n]];
    if(Rescale.Batches){
      original[[n]] = logcounts((nout[[n]][gene.chosen, ]))
    }else{
      original[[n]] = logcounts((sce[gene.chosen, which(sce$batches == bc.uniq[n])])) 
    }
    fscs = c(fscs, sce$FSC_log2[which(sce$batches == bc.uniq[n])])
  }
  
  # Slightly convoluted call to avoid re-writing code later.
  # Equivalent to fastMNN(GSE81076, GSE85241, k=20, d=50, approximate=TRUE)
  # Comments from Aaron: https://www.bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html#5_examining_the_effect_of_correction
  # The k= parameter specifies the number of nearest neighbours to consider when defining MNN pairs. This should be interpreted as the minimum frequency of each cell type or state in each batch.
  # Larger values will improve the precision of the correction by increasing the number of MNN pairs. It also provides some robustness to violations of the assumption that the batch vector is orthogonal to the biological subspace (Haghverdi et al. 2018), by allowing the neighbour search to ignore biological variation in each batch to identify the correct MNN pairs.
  # However, larger values of k can also reduce accuracy by allowing incorrect MNN pairs to form between cells of different types. Thus, we suggest starting with the default k and increasing it if one is confident that the same cell types are not adequately merged across batches. This is better than starting with a large k as incorrect merging is much harder to diagnose than insufficient merging.
  # When BSPARAM=IrlbaParam(deferred=TRUE), fastMNN() uses methods from the irlba package to perform the principal components analysis quickly. While the run-to-run differences should be minimal, it does mean that set.seed() is required to obtain fully reproducible results. The deferred= argument instructs fastMNN() to sacrifice some numerical precision for greater speed.
  set.seed(1001)
  mnn.out <- do.call(fastMNN, c(original, list(k=k.mnn, cos.norm = cos.norm, d=nb.pcs, auto.order=order2correct,
                                               approximate=TRUE)))
  dim(mnn.out$corrected)
  mnn.out$batch
  Rle(mnn.out$batch)
  #metadata(mnn.out)$merge.info$pairs[[1]]
  reducedDim(sce, "MNN") <- mnn.out$corrected
  sce$mnn_Batch <- as.character(mnn.out$batch)
  sce
  
  ##########################################
  # check the effectiveness of batch correction with MNN
  # 1) check the MNN pairs and lost variances
  # 2) visualization wiht PCA, and UMAP before and after correction
  # 3) kBET test
  ##########################################
  pdfname = paste0(resDir, "/scRNAseq_filtered_test_MNNbatchCorrection_effectiveness_v1.pdf")
  pdf(pdfname, width=14, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  sce <- runUMAP(sce, use_dimred="MNN", n_dimred = 15, ncomponents = 2)
  p = plotUMAP(sce, ncomponents = 2, colour_by="mnn_Batch", size_by = "FSC_log2", point_size= 0.01) + ggtitle("Corrected") 
  plot(p)
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
  p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") + 
    fontsize + ggtitle("Seurat clustering")
  p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
    fontsize + ggtitle("Seurat clustering")
  multiplot(p1, p2, cols = 2)
  
  sce = runTSNE(sce, use_dimred="MNN", perplexity = 20, n_dimred = 20, scale_features = TRUE)
  p1 = plotTSNE(sce, colour_by="pha-4", size_by = "FSC_log2") + 
    fontsize + ggtitle("Seurat clustering")
  p2 = plotTSNE(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
    fontsize + ggtitle("Seurat clustering")
  multiplot(p1, p2, cols = 2)
  
  # mnn.out$pairs
  source("scRNAseq_functions.R")
  Check.MNN.pairs(mnn.out, fscs)
  
  #omat <- do.call(cbind, original)
  #sce.qc <- SingleCellExperiment(list(logcounts=omat))
  set.seed(1000)
  with.var <- do.call(fastMNN, c(original,
                                 list(k=k.mnn, cos.norm = cos.norm, d=nb.pcs, auto.order=order2correct, approximate=TRUE,
                                      compute.variances=TRUE)))
  with.var$lost.var
  
  fsc.boundary = c(27.2, 28.5)
  bsc.boundary = c(24.5, 27)
  plotColData(sce,
              x = "FSC_log2",
              y = "BSC_log2",
              colour_by = "mnn_Batch",
              shape_by = "mnn_Batch"
              
  ) + geom_vline(xintercept = fsc.boundary, linetype="dashed", color = "blue", size=0.5) +
    geom_hline(yintercept= bsc.boundary, linetype="dashed", color = "red", size=0.5) 
      
  
  sel.tmp = which(sce$FSC_log2 > fsc.boundary[1] & sce$FSC_log2 < fsc.boundary[2] 
                  & sce$BSC_log2 > bsc.boundary[1] & sce$BSC_log2 < bsc.boundary[2])
  #sce.tmp = sce[gene.chosen, which(sce$mnn_Batch > 2)]
  sce.tmp = sce[, sel.tmp]
  #sce.tmp = sce
  dim(sce.tmp)
  
  sce.tmp <- runPCA(sce.tmp, ncomponents = 50, ntop=nrow(sce.tmp), method="irlba", exprs_values = "logcounts", scale_features = TRUE)
  
  plotPCA(sce.tmp, colour_by="batches") + ggtitle("Original")
  
  dff = as.data.frame(reducedDim(sce.tmp, "MNN"))
  colnames(dff) = paste0("PC", c(1:ncol(dff)))
  dff$batches = sce.tmp$batches
  ggp = ggplot(data=dff, aes(PC1, PC2, color=batches)) + geom_point(size=2) 
  plot(ggp);
  
  # Using irlba to set up the t-SNE, for speed.
  #set.seed(100)
  #osce <- runTSNE(sce.tmp, use_dimred="PCA", perplexity = 20)
  #ot <- plotTSNE(osce, colour_by="mnn_Batch") + ggtitle("Original")
  #set.seed(100)
  #csce <- runTSNE(sce, use_dimred="MNN", perplexity = 20)
  #ct <- plotTSNE(csce, colour_by="mnn_Batch") + ggtitle("Corrected")
  #multiplot(ot, ct, cols=2)
  
  set.seed(100)
  osce <- runUMAP(sce.tmp, use_dimred="PCA", perplexity = 20)
  ot <- plotUMAP(osce, colour_by="mnn_Batch", size_by = "FSC_log2") + ggtitle("Original")
  set.seed(100)
  csce <- runUMAP(sce.tmp, use_dimred="MNN", perplexity = 20)
  ct <- plotUMAP(csce, colour_by="mnn_Batch", size_by = "FSC_log2") + ggtitle("Corrected")
  multiplot(ot, ct, cols=2)
  
  kbet.orig <- kBET(
    df = t(reducedDim(sce.tmp, "PCA")), 
    batch = sce.tmp$mnn_Batch,
    heuristic = TRUE,
    do.pca = FALSE,
    verbose = TRUE, 
    addTest = FALSE,
    n_repeat = 200,
    plot = TRUE)
  
  #require('FNN')
  # data: a matrix (rows: samples, columns: features (genes))
  #k0=floor(mean(table(sce.tmp$batches))) #neighbourhood size: mean batch size 
  #knn <- get.knn(, k=k0, algorithm = 'cover_tree')
  kbet.bc = kBET(
    df = t(reducedDim(sce.tmp, "MNN")), 
    batch = sce.tmp$mnn_Batch,
    do.pca = FALSE,
    heuristic = FALSE,
    verbose = TRUE, 
    addTest = FALSE,
    n_repeat = 200,
    plot = TRUE)
  
  dev.off()
  
}

save(sce, file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_bcMNN.Rdata'))

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
library(scater)
library(scran)
library(scRNA.seq.funcs)
library(matrixStats)
#library(M3Drop)
library(RColorBrewer)
library(SingleCellExperiment)
#set.seed(100)
load(file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_bcMNN.Rdata')) 

Seurat.clustering = TRUE
Test.scran.clustering = FALSE
##########################################
# (test clustering method) but Seurat method is used at the end 
# https://master.bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/work-1-reads.html/
# DE analysis is tested with scran but will be considered to be replaced by MAST or Seurat
# DE analysis (or marker gene finding) following the cluster analysis

##########################################
sce <- runUMAP(sce, use_dimred="MNN", perplexity = 20, n_dimred = 20)
p = plotUMAP(sce, colour_by="mnn_Batch", size_by = "FSC_log2", point_size= 0.01) + ggtitle("Corrected") 
plot(p)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") + 
  fontsize + ggtitle("Seurat clustering")
p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
  fontsize + ggtitle("Seurat clustering")
multiplot(p1, p2, cols = 2)

if(Seurat.clustering)
{ 
  require(Seurat)
  library(cowplot)
  
  check.individualExample.geneMarker = FALSE
  ntops = 3 # nb of top gene markers
  resolutions = c(0.4, 0.6, 0.8, seq(1.0, 4.0, by = 0.5))
  
  for(rr in resolutions){
    pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_batchCorrected_clustering.Seurat_geneMarkers.scran_resolution", 
                     rr, ".pdf")
    pdf(pdfname, width=22, height = 18)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    ##########################################
    # test graph-based Louvain algorithm 
    ##########################################
    rr = 1.2
    
    pbmc = as.Seurat(sce)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    pbmc <- ScaleData(pbmc, features = rownames(pbmc))
    
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
    
    pbmc = FindNeighbors(object = pbmc, reduction = "MNN", k.param = 20, dims = 1:20)
    
    cat("--- resolution is :", rr, "---\n")
    pbmc = FindClusters(pbmc, resolution = rr, algorithm = 3)
    
    # DimPlot(pbmc, reduction = "umap")
    
    # note that you can set `label = TRUE` or use the LabelClusters function to help label
    # individual clusters
    pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10, reduction = "MNN", reduction.key = "umap", n.neighbors = 15, repulsion.strength = 2)
    DimPlot(pbmc, reduction = "umap")
    FeaturePlot(pbmc, features = c("pha-4", "hnd-1"), reduction = "umap")
    # FeaturePlot(pbmc, features = c("pha-4", "hnd-1"), reduction = "UMAP")
    
    sce$cluster_seurat <- factor(pbmc@active.ident)
    sce$cluster <- factor(pbmc@active.ident)
    sce = runPCA(sce)
    
    # configuration of umap
    # https://github.com/cran/umap/blob/master/R/umap.R
    umap.defaults = list(
      n_neighbors=20,
      n_components=2,
      metric="euclidean",
      n_epochs=200,
      input="data",
      init="spectral",
      min_dist=0.5,
      set_op_mix_ratio=1,
      local_connectivity=1,
      bandwidth=1.0,
      alpha=1,
      gamma=1.0,
      negative_sample_rate=5,
      a=NA,
      b=NA,
      spread=1,
      random_state=NA,
      transform_state=NA,
      knn_repeats=1,
      verbose=TRUE,
      umap_learn_args = NA
    )
    class(umap.defaults) = "umap.config"
    
    sce <- runUMAP(sce, use_dimred="MNN", perplexity = 20, n_dimred = 50, scale_features = TRUE, 
                   method = "umap-learn", config = umap.defaults)
    fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
    #plotUMAP(sce, colour_by="cluster", size_by = "FSC_log2") + 
    #  fontsize + ggtitle("Seurat clustering")
    p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") + 
      fontsize + ggtitle("Seurat clustering")
    p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
      fontsize + ggtitle("Seurat clustering")
    multiplot(p1, p2, cols = 2)
    
    
    sce = runTSNE(sce, use_dimred="MNN", perplexity = 20, n_dimred = 50, scale_features = FALSE)
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
                colour_columns_by=c("cluster"),
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

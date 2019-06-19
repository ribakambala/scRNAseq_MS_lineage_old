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
Norm.Vars.per.batch = TRUE
#rescaling.for.multBatch = TRUE
Use.fastMNN = TRUE
Test.mnnCorrect = FALSE

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

sce$FSC_log10 = log10(sce$FSC)
sce$BSC_log10 = log10(sce$BSC)
sce = sce[, which(sce$nb.cells == 1)]

plotColData(sce,
            x = "FSC_log10",
            y = "BSC_log10",
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
order2correct = c(3, 4, 2, 1)

if(!Norm.Vars.per.batch){
  ## here we use batches as blockes, i.e. calculated mean and variances separately for each batches and then fit the trend
  ## In doing so, we implicitly assume that the trend is the same between batches
  ## https://bioconductor.org/packages/3.10/workflows/vignettes/simpleSingleCell/inst/doc/var.html#41_using_the_block=_argument
  fit <- trendVar(sce, block=sce$batches, parametric=TRUE, assay.type="logcounts", use.spikes=FALSE)
  dec <- decomposeVar(sce, fit)
  
  par(mfrow=c(1,2))
  cols = c(1:length(bc.uniq)) 
  nb.genes = min(table(sce$batches))
  
  matplot(fit$means[c(1:nb.genes), ], fit$vars[c(1:nb.genes), ], col=cols,
          xlab="Mean log-expression", ylab="Variance of log-expression")
  curve(fit$trend(x), add=TRUE, col="red")
  
  plot(dec$mean, dec$total, xlab="Mean log-expression", 
       ylab="Variance of log-expression", pch=16)
  curve(fit$trend(x), col="dodgerblue", add=TRUE)
  
  tmp.sce <- sce
  tmp.sce$log_size_factor <- log(sizeFactors(sce))
  plotColData(tmp.sce, x="batches", y="log_size_factor")
  #p2 <- plotColData(tmp.416B, x="Plate", y="log_size_factor_ERCC")
  #multiplot(p1, p2, cols=2)
  
  dec.sorted <- dec[order(dec$bio, decreasing=TRUE),]
  head(dec.sorted)
  
  # here HVGs selected with FDR<0.01
  gene.chosen <- which(dec.sorted$FDR < 0.05)
  length(gene.chosen)

}else{
  ## here we search batch-specific features
  ## https://bioconductor.org/packages/3.10/workflows/vignettes/simpleSingleCell/inst/doc/var.html#42_fitting_batch-specific_trends
  sce.bc = sce
  ## recommendation in the function help: Centering should be performed by running multiBlockNorm before calling this function. 
  sce.bc <- multiBlockNorm(sce.bc, block = sce.bc$batches)
  
  par(mfrow=c(1,1))
  plot(sizeFactors(sce), sizeFactors(sce.bc), log='xy'); abline(0, 1, lwd =2, col = 'red') # did not change anything here, weired
  
  comb.out <- multiBlockVar(sce.bc, block=sce.bc$batches, assay.type="logcounts", trend.args=list(parametric=TRUE, use.spikes=FALSE))
  
  #comb.out = multiBlockVar(sce.bc, block = sce.bc$batches, trend.args = list(use.spikes = FALSE))
  head(comb.out[,1:6])
  
  par(mfrow=c(1, length(bc.uniq)))
  #is.spike <- isSpike(sce.416B.2)
  for (plate in unique(sce.bc$batches)) {
    cur.out <- comb.out$per.block[[plate]]
    plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
         ylab="Variance of log-expression", main=plate, ylim = c(0, 25), xlim = c(0, 15))
    curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
    #points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
  }
  
  dec.bc.sorted <- comb.out[order(comb.out$bio, decreasing=TRUE), ]
  head(dec.bc.sorted)
  
  #length(which(dec.sorted$bio>0))
  # here HVGs selected with FDR<0.01
  #gene.chosen.bc <- which(dec.bc.sorted$p.value < 0.05)
  #gene.chosen.bc <- which(dec.bc.sorted$FDR < 0.1)
  gene.chosen.bc = which(dec.bc.sorted$bio>0)
  #length(which(dec.sorted$bio>0))
  length(gene.chosen.bc)
  
  ### compare the block and batch-specific HVGs selection
  #cat(length(gene.chosen), length(gene.chosen.bc), length(intersect(gene.chosen, gene.chosen.bc)), "\n")
  #library(VennDiagram)
  #venn.diagram(
  #  x = list(gene.chosen, gene.chosen.bc),
  #  category.names = c("batch-block" , "batch-specific"), filename = paste0(resDir, "/batch_block_specific.png"), output = TRUE)
  
  gene.chosen = gene.chosen.bc
  
}
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
  
  Rescale.Batches = FALSE
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
    fscs = c(fscs, sce$FSC_log10[which(sce$batches == bc.uniq[n])])
  }
  
  # Slightly convoluted call to avoid re-writing code later.
  # Equivalent to fastMNN(GSE81076, GSE85241, k=20, d=50, approximate=TRUE)
  # Comments from Aaron: https://www.bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html#5_examining_the_effect_of_correction
  # The k= parameter specifies the number of nearest neighbours to consider when defining MNN pairs. This should be interpreted as the minimum frequency of each cell type or state in each batch.
  # Larger values will improve the precision of the correction by increasing the number of MNN pairs. It also provides some robustness to violations of the assumption that the batch vector is orthogonal to the biological subspace (Haghverdi et al. 2018), by allowing the neighbour search to ignore biological variation in each batch to identify the correct MNN pairs.
  # However, larger values of k can also reduce accuracy by allowing incorrect MNN pairs to form between cells of different types. Thus, we suggest starting with the default k and increasing it if one is confident that the same cell types are not adequately merged across batches. This is better than starting with a large k as incorrect merging is much harder to diagnose than insufficient merging.
  # When BSPARAM=IrlbaParam(deferred=TRUE), fastMNN() uses methods from the irlba package to perform the principal components analysis quickly. While the run-to-run differences should be minimal, it does mean that set.seed() is required to obtain fully reproducible results. The deferred= argument instructs fastMNN() to sacrifice some numerical precision for greater speed.
  set.seed(1001)
  mnn.out <- do.call(fastMNN, c(original, list(k=15, cos.norm = TRUE, d=50, auto.order=order2correct,
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
  pdfname = paste0(resDir, "/scRNAseq_filtered_test_MNNbatchCorrection_effectiveness.pdf")
  pdf(pdfname, width=18, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  # mnn.out$pairs
  source("scRNAseq_functions.R")
  Check.MNN.pairs(mnn.out, fscs)
  
  sce <- runUMAP(sce, use_dimred="MNN", perplexity = 20)
  #sce$FSC_log10 = sce$FSC_log10 - 5
  plotUMAP(sce, colour_by="mnn_Batch", size_by = "FSC_log10", point_size= 0.01) + ggtitle("Corrected") 
  
  #omat <- do.call(cbind, original)
  #sce.qc <- SingleCellExperiment(list(logcounts=omat))
  set.seed(1000)
  with.var <- do.call(fastMNN, c(original,
                                 list(k=20, cos.norm = TRUE, d=50, auto.order=order2correct, approximate=TRUE,
                                      compute.variances=TRUE)))
  with.var$lost.var
  
  plotColData(sce,
              x = "FSC_log10",
              y = "BSC_log10",
              colour_by = "mnn_Batch",
              shape_by = "mnn_Batch"
              
  ) + geom_hline(yintercept= c(4.9, 5.3)  , linetype="dashed", color = "darkgray", size=0.5) +
      geom_vline(xintercept = c(5.45, 5.75), linetype="dashed", color = "black", size=0.5)
  
  sel.tmp = which(sce$BSC_log10>4.9 & sce$BSC_log10 < 5.3 & sce$FSC_log10>5.45 & sce$FSC_log10<5.75)
  #sce.tmp = sce[gene.chosen, which(sce$mnn_Batch > 2)]
  sce.tmp = sce[, sel.tmp]
  sce.tmp = sce
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
  ot <- plotUMAP(osce, colour_by="mnn_Batch", size_by = "FSC_log10") + ggtitle("Original")
  set.seed(100)
  csce <- runUMAP(sce.tmp, use_dimred="MNN", perplexity = 20)
  ct <- plotUMAP(csce, colour_by="mnn_Batch", size_by = "FSC_log10") + ggtitle("Corrected")
  multiplot(ot, ct, cols=2)
  
  kbet.orig <- kBET(
    df = t(reducedDim(sce.tmp, "PCA")), 
    batch = sce.tmp$batches,
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

if(Use.mnnCorrect){
  set.seed(100)
  require(BiocParallel)
  require(stats)
  bpp <- MulticoreParam(5)
  bpp
  
  mnn.out2 = mnnCorrect(R6879, R7116, R7130, k = 20, sigma = 0.1, cos.norm.in = TRUE, cos.norm.out = TRUE, order = c(3, 2, 1), 
                        svd.dim = 0, var.adj = FALSE, subset.row = gene.chosen, pc.approx = TRUE, BPPARAM=bpp)
  
  head(mnn.out2$pairs[[2]])
  
  ## check the pairs of cells used between batches
  Check.SNN.pairs = FALSE
  if(Check.SNN.pairs){
    for(n in 2:length(mnn.out$pairs))
    {
      #n = 3
      paires = data.frame(mnn.out2$pairs[[n]])
      
      cat("cell in the batch", n, ":  ",
          "total pairs --", nrow(paires), ",",  
          "paired cell in batch", n-1, "-- ", length(unique(paires[which(paires[,3]==(n-1)), 2])), ",", 
          "paired cell in batch", n+1, "-- ", length(unique(paires[which(paires[,3]==(n+1)), 2])), "\n"
      )
    }
  }
  
  res1 <- mnn.out2$corrected[[1]]
  res2 <- mnn.out2$corrected[[2]]
  res3 <- mnn.out2$corrected[[3]]
  
  dimnames(res1) <- dimnames(R6879)
  dimnames(res2) <- dimnames(R7116)
  dimnames(res3) <- dimnames(R7130)
  res = cbind(res1, res2, res3)
  
  assay(sce, "corrected") <- res
  
  #sce = sce.qc
  #save(sce, mnn.out, mnn.out2, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_batchCorrectMNN_SCE.Rdata'))
  
  ##########################################
  # double check the relatinship between fastMNN and mnnCorrect 
  ##########################################
  #osce = runPCA(sce, ncomponents = 50, ntop=Inf, method="irlba", exprs_values = "logcounts", feature_set = gene.chosen)
  csce <- runPCA(sce, ncomponents = 50, method="irlba", exprs_values = "corrected", feature_set = gene.chosen,
                 scale_features = TRUE, detect_outliers = FALSE)
  
  xx = as.data.frame(reducedDim(csce, "MNN"));
  yy = as.data.frame(reducedDim(csce, "PCA"))
  
  head(xx[, c(1:10)])
  head(yy[, c(1:10)])
  
  par(mfrow = c(1, 2))
  plot(xx[, c(1:2)], xlab = 'PC1', ylab = "PC2", main = "lowDim output from fastMNN")
  plot(yy[, c(1:2)], xlab = 'PC1', ylab = "PC2", main = "PCA from mnnCorret output")
  #plot(xx[, c(3:4)], xlab = 'PC1', ylab = "PC2", main = "lowDim output from fastMNN")
  #plot(yy[, c(3:4)], xlab = 'PC1', ylab = "PC2", main = "PCA from mnnCorret output")
  
  
  pdfname = paste0(resDir, "/scRNAseq_filtered_test_batchCorrection_fastMNNlowDim_vs_mnnCorrectOutput.pdf")
  pdf(pdfname, width=18, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  od = scater::plotPCA(
    osce,
    run_args = list(exprs_values = "logcounts"), 
    size_by = "total_counts",
    #size_by = "total_features_by_counts",
    colour_by = "Batch"
  ) + ggtitle(paste0("PCA -- origine "))
  
  cd = scater::plotPCA(
    csce,
    run_args = list(exprs_values = "corrected"),
    size_by = "total_counts",
    #size_by = "total_features_by_counts",
    colour_by = "request"
  ) + ggtitle(paste0("PCA -- corrected by mnn"))
  
  multiplot(od, cd, cols=2)
  
  dev.off()
  
}

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
load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_batchCorrectMNN_SCE.Rdata'))

library(scater)
library(scran)
library(scRNA.seq.funcs)
library(matrixStats)
#library(M3Drop)
library(RColorBrewer)
library(SingleCellExperiment)
#set.seed(100)

##########################################
# test clustering methods in scran 
# https://master.bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/work-1-reads.html/
##########################################
TEST.Aaron.scRNAseq.workflow.clustering.part = FALSE
if(TEST.Aaron.workflow)
{ 
  sce = runPCA(sce, ncomponents = 50, ntop=Inf, method="irlba", exprs_values = "corrected")
  
  #set.seed(100)
  #sce <- runTSNE(sce, use_dimred="MNN", perplexity = 20, n_dimred = 20)
  
  set.seed(100)
  sce <- runTSNE(sce, use_dimred="MNN", perplexity = 20, n_dimred = 20)
  
  set.seed(100)
  sce = runUMAP(sce, use_dimred="MNN", perplexity = 20, n_dimred = 20)
  
  #sce = runDiffusionMap(sce, use_dimred = "MNN", n_dimred = 20)
  
  pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_batchCorrected_clustering_testing.pdf")
  pdf(pdfname, width=10, height = 6)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
  
  ##########################################
  # test grapha method from scran
  ##########################################
  snn.gr <- buildSNNGraph(sce, use.dimred="MNN")
  clusters <- igraph::cluster_walktrap(snn.gr)
  table(clusters$membership, sce$Batch)
  
  sce$cluster <- factor(clusters$membership)
  plotTSNE(sce, colour_by="cluster") + ggtitle("scran -- graph based clustering")
  plotUMAP(sce, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
    fontsize + ggtitle("scran -- graph based clustering")
  
  ##########################################
  # test hierachy clustering from scran
  ##########################################
  pcs <- reducedDim(sce, "MNN")
  my.dist <- dist(pcs)
  my.tree <- hclust(my.dist, method="ward.D2")
  
  #hist(my.tree)
  library(dynamicTreeCut)
  my.clusters <- unname(cutreeDynamic(my.tree, cutHeight = 3, 
                                      distM=as.matrix(my.dist), 
                                      minClusterSize=5, verbose=0))
  
  table(my.clusters, sce$Batch)
  
  sce$cluster <- factor(my.clusters)
  
  plotTSNE(sce, colour_by="cluster", size_by = "total_features_by_counts") + fontsize + ggtitle("scran -- hcluster")
  plotUMAP(sce, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
    fontsize + ggtitle("scran -- hcluster")
  
  #plotDiffusionMap(sce, colour_by="cluster", size_by = "total_features_by_counts") + fontsize
  
  library(cluster)
  clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
  sil <- silhouette(my.clusters, dist = my.dist)
  sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
  sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
  plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
       border=sil.cols, col=sil.cols, do.col.sort=FALSE) 
  
  
  ##########################################
  # test graph-based Louvain algorithm 
  ##########################################
  require(Seurat)
  library(cowplot)
  #srt = Seurat::Convert(from = sce, to = "seurat") 
  pbmc = as.Seurat(sce)
  
  #Seurat::DimPlot(pbmc, dims = c(1, 2), reduction = "MNN")
  pbmc = FindNeighbors(object = pbmc, reduction = "MNN", k.param = 10, dims = 1:10)
  pbmc = FindClusters(pbmc, resolution = 1, algorithm = 3)
  sce$cluster_seurat <- factor(pbmc@active.ident)
  sce$cluster <- factor(pbmc@active.ident)
  
  plotTSNE(sce, colour_by="cluster", size_by = "total_features_by_counts") + fontsize + ggtitle("seurat - graph base clustering")
  plotUMAP(sce, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
    fontsize + ggtitle("seurat -- graph based clustering")
  
  dev.off()
  
}
##########################################
# DE analysis (or marker gene finding) following the cluster analysis
# To do it, we also have several options
# 1) scran 
##########################################
Find.Gene.Markers.with.scran = FALSE
if(Find.Gene.Markers.with.scran){
  my.clusters = as.numeric(as.character(sce$cluster_seurat))
  
  design <- model.matrix( ~ sce$Batch)
  design <- design[,-1,drop=FALSE]
  
  # run the find markers and then collect markers for each clusters
  markers <- findMarkers(sce, my.clusters, design = design)
  
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
  
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
  
  pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_batchCorrected_clustering_markerGenes_examples.pdf")
  pdf(pdfname, width=16, height = 12)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  plotTSNE(sce, colour_by="cluster", size_by = "total_features_by_counts") + 
    ggtitle("seurat - graph base clustering") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
  
  plotUMAP(sce, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
    fontsize + ggtitle("seurat -- graph based clustering")
  
  plotHeatmap(sce, features=top.markers,
              columns=order(sce$cluster_seurat), 
              colour_columns_by=c("cluster"),
              cluster_cols=FALSE, show_colnames = FALSE,
              center=TRUE, symmetric=TRUE, zlim=c(-5, 5))
  
  for(n in 1:length(top.markers)) {
    xx = plotTSNE(sce, colour_by = top.markers[n]) 
    plot(xx)
  }
  
  dev.off()
  
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

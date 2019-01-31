########################################################
########################################################
# Section : Trojectory Analysis
# The methods to test are following:
# 1) Monocle2
# 2) Slingshot
# 3) TSCAN
# 4) Celltree
# 5) PAGA
# 5) ELPiGrapha (possible)
# destiny is installed but not to test in the first time
# imputation method MAGIC is also installed
########################################################
########################################################
Test.Monocole = FALSE
if(Test.Monocole){
  library(monocle)
  library(dplyr)
  library(reshape2)
  
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata')) 
  
  # importCDS(sce, import_all = FALSE)
  HSMM <- newCellDataSet(as.matrix(counts(sce)),
                         expressionFamily=VGAM::negbinomial.size())
  
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  print(head(fData(HSMM)))
  expressed_genes <- row.names(subset(fData(HSMM),
                                      num_cells_expressed >= 10))
  
  L <- log(exprs(HSMM[expressed_genes,]))
  melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
  
  qplot(value, geom = "density", data = melted_dens_df) +
    stat_function(fun = dnorm, size = 0.5, color = 'red') +
    xlab("Standardized log(FPKM)") +
    ylab("Density")
  
  disp_table <- dispersionTable(HSMM)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
  plot_ordering_genes(HSMM)
  
  plot_pc_variance_explained(HSMM, return_all = F, max_components = 40, norm_method = 'log') # norm_method='log'
  
  HSMM <- reduceDimension(HSMM, max_components = 2, norm_method = c("log"), 
                          reduction_method = 'ICA', verbose = T)
  
  HSMM <- clusterCells(HSMM, method = "DDRTree")
  plot_cell_clusters(HSMM, 1, 2, color = "Cluster")
  
  
  HSMM_myo = HSMM
  diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,], 
                                        fullModelFormulaStr = "~ Cluster", 
                                        reducedModelFormulaStr = "~1" 
  )
  ordering_genes <- row.names (diff_test_res[order(diff_test_res$pval), ])[1:200]
  
  HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
  plot_ordering_genes(HSMM_myo)
  
  HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,
                              method = 'DDRTree')
  
  HSMM_myo <- orderCells(HSMM_myo)
  plot_cell_trajectory(HSMM_myo, color_by = "Cluster")
  
}

Test.cellTree = FALSE
if(Test.cellTree){
  library(cellTree)
  library(scran)
  library(scater)
  #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata')) 
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  
  lda.results = compute.lda(as.matrix(logcounts(sce.HVG.Brenneck)), k.topics=3:8, method="maptpx", log.scale = FALSE)
  # lda.results = compute.lda(as.matrix(logcounts(sce)), k.topics=6, method="Gibbs", log.scale = FALSE)
  
  mst.tree = compute.backbone.tree(HSMM_lda_model, only.mst=TRUE, rooting.method = "longest.path", width.scale.factor = 2.0)
  # plot the tree (showing topic distribution for each cell):
  mst.tree.with.layout = ct.plot.topics(mst.tree)
  
}

Test.Slingshot = FALSE
if(Test.Slingshot){
  #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata')) 
  #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  
  
  sim = sce.HVG.Brenneck
  
  library(slingshot, quietly = FALSE)
  
  assay(sim, "norm") <- exp(logcounts(sim))
  
  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }
  
  #assays(sim)$norm <- FQnorm(assays(sim)$counts)
  
  pca <- prcomp(t((assays(sim)$logcounts)), scale. = FALSE, center = TRUE)
  #xx = t((assays(sce)$logcounts))
  #vars = apply(xx, 2, sd)
  #ntop = 1000; 
  #xx = xx[, order(vars, decreasing = TRUE)]
  #xx = xx[, c(1:ntop)]
  #pca = prcomp(xx, scale. = TRUE, center = TRUE)
  
  rd1 <- pca$x[,1:2]
  plot(rd1, col = 'red', pch=16)
  
  
  library(destiny, quietly = TRUE)
  dm <- DiffusionMap(t(log(assays(sim)$norm)))
  rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)
  
  plot(rd2, col = topo.colors(100), pch=16, asp = 1)
  
  reducedDims(sim) <- SimpleList(PCA = rd1, DiffMap = rd2)
  
  library(mclust, quietly = TRUE)
  
  cl1 <- Mclust(rd2)$classification
  colData(sim)$GMM <- cl1
  
  library(RColorBrewer)
  plot(rd2, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
  
  
  sce <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'DiffMap')
  
  summary(sce$slingPseudotime_1)
  
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  plot(reducedDims(sce)$DiffMap, col = colors[cut(sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2)
  
}


Test.ELPiGraph = FALSE
if(Test.ELPiGraph){
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  
  library(ElPiGraph.R)
  TreeEPG <- computeElasticPrincipalTree(X = t(logcounts(sce.HVG.Brenneck)), NumNodes = 7, NumEdges=6, Lambda = .03, Mu = .01)
  
}


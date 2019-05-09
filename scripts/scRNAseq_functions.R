##################################################
##################################################
## Project: collections of functions for single-cell RNAseq data analysis
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Thu Feb 22 14:51:03 2018
##################################################
##################################################
process.countTable = function(all, design, additional.filter = NULL)
{
  newall = data.frame(as.character(all[,1]), stringsAsFactors = FALSE)
  
  for(n in 1:nrow(design))
  {
    #n = 1;
    ## found the matching ID in design matrix
    if(!is.null(additional.filter)){
      jj = intersect(grep(design$SampleID[n], colnames(all)), grep(additional.filter, colnames(all)));
    }else{
      jj = grep(design$SampleID[n], colnames(all));
    }
    
    ## deal with several mapping 
    if(length(jj)==1) {
      #index = c(index,jj)
      newall = data.frame(newall, all[, jj])
    }else{
      cat(length(jj), " samples found for ID", design$SampleID[n], "\n")
      cat("start to merge those samples considering them as technical replicates...\n")
      newall = data.frame(newall, apply(as.matrix(all[, jj]), 1, sum))
    }
  }
  
  colnames(newall)[1] = "gene";
  jj = which(colnames(design)=="SampleID")
  o1 = c(setdiff(c(1:ncol(design)), jj),jj)
  #design = design
  colnames(newall)[-1] = apply(design[, o1], 1, function(x) paste0(x, collapse = "_"))
  
  #if(time.series){
  #  colnames(newall)[-1] = paste0(design$stage, "_", design$treatment, "_", design$SampleID)
  #}else{
  #}
  
  return(newall)
}

aggrate.nf.QCs = function(dirs.all, modules2cat = c("star", "featureCounts"))
{
  # dir = '../../../Ariane/R7116_R7130_scrnaseq/results_all/multiqc_data_1'
  # modules2cat = c("star", "featureCounts")
  
  res = NULL;
  for(fold in dirs.all){
    if(dir.exists(fold)){
      
      yy = NULL;
      for (mod in modules2cat)
      {
        # fold = dirs.all[1]; mod = modules2cat[1];
        files = list.files(path = fold, pattern = "*.txt", full.names = TRUE)
        files = files[grep(mod, files)]
        if(length(files)==0){
          cat("Error -- no file found for: ",  mod, "in dir --", fold, "\n")
        }else{
          for(f in files){
            if(is.null(yy)){
              yy = data.frame(read.delim(f, sep = "\t", header = TRUE), stringsAsFactors = FALSE)
            }else{
              test = data.frame(read.delim(f, sep = "\t", header = TRUE), stringsAsFactors = FALSE)
              yy = merge(yy, test, by = "Sample", suffixes = c("",".1"))
            }
          } 
        }
      } 
      
    }else{
      cat('Error -- ', fold, "not exist--\n")
    }
    
    if(is.null(res)){
      res = yy;
    }else{
      res = rbind(res, yy)
    }
  }
  
  return(res)
}


##########################################
# function for technical replicates 
##########################################
compare.techinical.replicates = function(design, counts)
{
  library(SingleCellExperiment)
  library(scater)
  options(stringsAsFactors = FALSE)
  
  # change seqInfos in design to label technicalreps 
  design$seqInfos[which(design$seqInfos=="R7130_HLWTCBGX9_1")] = "R7130_HHG5KBGX9_1_techrep_nexseq"
  design$seqInfos[which(design$seqInfos=="R7130_CCYTEANXX_4")] = "R7130_HHGHNBGX9_1_techrep_hiseq"
  design$seqInfos[which(design$seqInfos=="R7133_CD2GTANXX_5")] = "R7130_HHGHNBGX9_1_techrep_hiseq_R7133"
  
  ## add some new features for design for quality controls
  design$log10_Total = log10(design$total_reads)
  #design$percent_mapped = design$uniquely_mapped/design$T
  #design$percent_assigned = design$Assigned/design$uniquely_mapped
  design$percent_rRNA = design$rRNA / design$Total
  
  ## make SCE object and remove genes with zero reads detected
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = as.data.frame(design), 
                              rowData = data.frame(gene_names = rownames(counts), feature_symbol = rownames(counts)))
  #write.csv(counts(sce), file=paste0(tabDir, "scRNAseq_raw_readCounts", version.analysis, ".csv"), row.names=TRUE)
  keep_feature <- rowSums(counts(sce) > 0) > 0
  sce <- sce[keep_feature, ]
  
  
  #is.spike <- grepl("^ERCC", rownames(sce))
  is.mito <- rownames(sce) %in% gg.Mt;
  is.ribo <- rownames(sce) %in% gg.ribo;
  summary(is.mito)
  summary(is.ribo)
  
  sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito, Ribo=is.ribo))
  
  head(colnames(colData(sce)), 20)
  
  # some general statistics for each request and lane
  #cols = c(rep('gray', 2), 'red', 'darkblue', 'darkred', 'blue', 'red')
  plotColData(sce, y = "log10_Total", x = "seqInfos") + ggtitle("total nb of reads")
  plotColData(sce, y="uniquely_mapped_percent", x="seqInfos") + ggtitle("% of uniquely mapped ")
  plotColData(sce, y="percent_assigned", x="seqInfos") + ggtitle("% of assigned")
  plotColData(sce, y="pct_counts_Ribo", x="seqInfos") + ggtitle("% of rRNA contamination")
  plotColData(sce, y="pct_counts_Mt", x="seqInfos") + ggtitle("% of Mt")
  
  plotColData(sce, y="log10_total_counts", x="seqInfos") + ggtitle("total nb of reads mapped to transcripts")
  plotColData(sce, y="total_features_by_counts", x="seqInfos") + ggtitle("total nb of genes")
  
  plotColData(sce, 
              x = "log10_Total",
              y = "uniquely_mapped_percent",
              #colour_by = "uniquely_mapped_percent",
              colour_by = "seqInfos",
              size_by = "pct_counts_Ribo"
  )
  
  plotColData(sce,
              x = "log10_total_counts",
              y = "log10_total_features_by_counts",
              #colour_by = "percent_mapped",
              colour_by = "seqInfos",
              size_by = "pct_counts_Mt"
  ) + scale_x_continuous(limits=c(4, 7)) +
    scale_y_continuous(limits = c(2.5, 4.1)) +
    geom_hline(yintercept=log10(c(500, 1000, 5000)) , linetype="dashed", color = "darkgray", size=0.5) +
    geom_vline(xintercept = c(4:6), linetype="dotted", color = "black", size=0.5)
  
  # filter cell and genes
  threshod.total.counts.per.cell = 10^4
  threshod.nb.detected.genes.per.cell = 1000;
  
  filter_by_total_counts <- (sce$total_counts > threshod.total.counts.per.cell)
  table(filter_by_total_counts)
  filter_by_expr_features <- (sce$total_features_by_counts > threshod.nb.detected.genes.per.cell)
  table(filter_by_expr_features)
  filter_by_MT = sce$pct_counts_Mt < 7.5
  table(filter_by_MT)
  
  sce$use <- (
    filter_by_expr_features & # sufficient features (genes)
      filter_by_total_counts & # sufficient molecules counted
      # filter_by_ERCC & # sufficient endogenous RNA
      filter_by_MT # remove cells with unusual number of reads in MT genes
  )
  table(sce$use)
  sce = sce[, sce$use]
  
  num.cells <- nexprs(sce, byrow=TRUE)
  ave.counts <- calcAverage(sce)
  genes.to.keep <- num.cells > 5 & ave.counts >= 1  & ave.counts <10^6  # detected in >= 2 cells, ave.counts >=5 but not too high
  summary(genes.to.keep)
  # remove mt and ribo genes
  genes.to.keep = genes.to.keep & ! rownames(sce) %in% gg.Mt & ! rownames(sce) %in% gg.ribo
  summary(genes.to.keep)
  
  sce <- sce[genes.to.keep, ]
  
  
  # select cells having technical replicates and normalize them  
  library(scRNA.seq.funcs)
  library(scater)
  library(scran)
  options(stringsAsFactors = FALSE)
  
  sce$sels = (sce$request !="R6875" & sce$request != "R7116" & sce$seqInfos != "R7130_HHG5KBGX9_1_techrep_nexseq" &
                sce$seqInfos != "R7130_HHG5KBGX9_1")
  sce.qc = sce[, sce$sels]
  
  reducedDim(sce.qc) <- NULL
  endog_genes <- !rowData(sce.qc)$is_feature_control
  
  set.seed(1234567)
  assay(sce.qc, "logcounts") <- log2(calculateCPM(sce.qc, use_size_factors = FALSE) + 1)
  #sizeFactors(sce.qc) = calculate.sizeFactors.DESeq2(counts(sce.qc))
  #sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
  #qclust <- quickCluster(sce.qc, min.size = 30)
  #sce.qc <- computeSumFactors(sce.qc, sizes = 15, clusters = qclust)
  #sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
  
  main = "cpm"
  scater::plotPCA(
    sce.qc[endog_genes, ],
    run_args = list(exprs_values = "logcounts"), 
    size_by = "total_counts",
    #size_by = "total_features_by_counts",
    colour_by = "seqInfos"
  ) + ggtitle(paste0("PCA -- ", main))
  
  set.seed(1234567)
  param.perplexity = 10;
  plotTSNE(
    sce.qc[endog_genes, ],
    run_args = list(exprs_values = "logcounts", perplexity = param.perplexity), 
    size_by = "total_counts",
    #size_by = "total_features_by_counts",
    colour_by = "seqInfos"  
  ) + ggtitle(paste0("tSNE - perplexity = ", param.perplexity, "--", main))
  
  plotUMAP(
    sce.qc[endog_genes, ],
    run_args = list(exprs_values = "logcounts"), 
    size_by = "total_counts",
    #size_by = "total_features_by_counts",
    colour_by = "seqInfos"
  ) + ggtitle(paste0("UMAP -- ", main))
  
  ## check the correction of the same cells from different technical replicates
  bcs = unique(sce.qc$barcodes)
  correlations = c()
  for(n in 1:length(bcs))
  {
    #n = 1
    xx = as.data.frame(logcounts(sce.qc[, which(sce.qc$barcodes == bcs[n])]))
    if(ncol(xx) == 3) correlations = rbind(correlations, c(cor(xx[, 1], xx[, 2]), cor(xx[, 1], xx[, 3]), cor(xx[, 2], xx[, 3])))
  }
  
  colnames(correlations) = c('rep0.vs.hiseq.rep1', 'rep0.vs.hiseq.rep2', 'hiseq.rep1.vs.hiseq.rep_2')
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

panel.fitting = function (x, y, bg = NA, pch = par("pch"), cex = 0.2, col='black') 
{
  #x = yy[,1];y=yy[,2];
  #kk = which(x>0 & y>0); x=x[kk];y=y[kk]
  lims = range(c(x, y), na.rm = TRUE)
  points(x, y, pch = 1, col = col, cex = cex, xlim=lims, ylim=lims)
  abline(0, 1, lwd=1.5, col='red')
  R = cor(x, y, use="na.or.complete", method='pearson')
  text(lims[2]*0.2, lims[2]*0.9, paste0('R = ', signif(R, d=2)), cex=1., col='red')
  #jj = which(!is.na(x) & !is.na(y))
  #fit = lm(y[jj] ~ x[jj])
  #slope=summary(fit)$coefficients[1]
  #slope = fit$coefficients[2]
  #intercept = fit$coefficients[1]
  #pval=summary(fit)$coefficients[4]
  #abline(intercept, slope, lwd=1.2, col='darkblue', lty=3)
  #text(lims[2]*0.1, lims[2]*0.7, paste0('slop = ', signif(slope, d=2)), cex=1., col='blue')
  #text(lims[2]*0.1, lims[2]*0.6, paste0('pval = ', signif(pval, d=2)), cex=1., col='blue')
  #ok <- is.finite(x) & is.finite(y)
  #if (any(ok)) 
  #lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
  #        col = col.smooth, ...)
}

plot.pair.comparison.plot = function(xx, linear.scale = TRUE, main = "pairwise.comparision"){
  yy = as.matrix(xx)
  if(linear.scale){
    yy[which(yy==0)] = NA;
    yy = log2(yy)
  }
  pairs(yy, lower.panel=NULL, upper.panel=panel.fitting, main = main)
  
}


##########################################
# several common functions for normalizations 
# from Hemberg lab
# hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#normalization-theory
##########################################
calc_cpm <- function (expr_mat, spikes = NULL) 
{
  norm_factor <- colSums(expr_mat[-spikes, ])
  return(t(t(expr_mat)/norm_factor)) * 10^6
}

Down_Sample_Matrix <- function (expr_mat)
{
  min_lib_size <- min(colSums(expr_mat))
  down_sample <- function(x) {
    prob <- min_lib_size/sum(x)
    return(unlist(lapply(x, function(y) {
      rbinom(1, y, prob)
    })))
  }
  down_sampled_mat <- apply(expr_mat, 2, down_sample)
  return(down_sampled_mat)
}

cal_uq_Hemberg = function (expr_mat, spikes = NULL) 
{
  UQ <- function(x) {
    quantile(x[x > 0], 0.75)
  }
  if(!is.null(spikes)){
    uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
  }else{
    uq <- unlist(apply(expr_mat, 2, UQ))
  }
  
  norm_factor <- uq/median(uq)
  return(t(t(expr_mat)/norm_factor))
}

calculate.sizeFactors.DESeq2 = function(expr_mat)
{
  # expr_mat = counts(sce.qc)
  require('DESeq2')
  condition <- factor(rep("A", ncol(expr_mat)))
  dds <- DESeqDataSetFromMatrix(expr_mat, DataFrame(condition), design = ~ 1)
  dds <- estimateSizeFactors(dds)
  
  return(sizeFactors(dds))
  
}





########################################################
########################################################
# Section :
# test code and not used anymore 
########################################################
########################################################
##########################################
# test seurat for clustering 
# original codes in https://satijalab.org/seurat/pbmc3k_tutorial.html
##########################################
Test.Seurat.workflow = FALSE
if(Test.Seurat.workflow){
  library(Seurat)
  library(dplyr)
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))
  
  #pbmc = as.seurat(from = sce)
  #pbmc = Convert(from = sce, to = 'seurat', raw.data.slot = "counts",data.slot = "logcounts")
  pbmc <- CreateSeuratObject(raw.data = counts(sce), min.cells = 0, min.genes = 200, 
                             project = "seurat_test")
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
  
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.5, y.cutoff = 0.5)
  
  length(x = pbmc@var.genes)
  
  pbmc <- ScaleData(object = pbmc)
  
  pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
                 genes.print = 5)
  
  PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
  VizPCA(object = pbmc, pcs.use = 1:2)
  
  PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
  
  pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
  
  PCHeatmap(object = pbmc, pc.use = 1:4, cells.use = 80, do.balanced = TRUE, label.columns = FALSE)
  
  pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
  JackStrawPlot(object = pbmc, PCs = 1:12)
  
  PCElbowPlot(object = pbmc)
  
  pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                       resolution = 2, print.output = 0, save.SNN = TRUE)
  
  PrintFindClustersParams(object = pbmc)
  #PrintCalcParams(pbmc, calculation = "FindClusters")
  
  pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE, perplexity=10, eta=2000)
  TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 3.0)
  
  ## test phate plot
  pbmc_phate <- RunPHATE(object = pbmc, gamma=1, npca = 20, k=5, seed.use = 10, plot.optimal.t = TRUE, t=14)
  # Plot results
  DimPlot(object = pbmc_phate, reduction.use = 'phate', pt.size = 2.0)
  #DimPlot(object = pbmc_phate, reduction.use = 'pca')
}


##################################################
# test Hamberg's single-cell RNA seq analysis, clustering part
# original code https://hemberg-lab.github.io/scRNA.seq.course/biological-analysis.html
##################################################
### test SC3
TEST.Hamberg.workflow.clustering = FALSE
if(TEST.Hamberg.workflow.clustering)
{
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata')) 
  library(SC3)
  #library()
  expr_matrix =  exp(logcounts(sce))
  Brennecke_HVG <- BrenneckeGetVariableGenes(
    expr_mat = expr_matrix,
    spikes = NA,
    fdr = 0.2,
    minBiolDisp = 0.2
  )
  #HVG_genes <- Brennecke_HVG
  
  ## another method to identify by Kolodziejczyk AA, Kim JK, Tsang JCH et al. (2015)
  assay(sce, "normcounts") <- exp(logcounts(sce))
  means <- rowMeans(normcounts(sce))
  cv2 <- apply(normcounts(sce), 1, var)/means^2
  dm.stat <- DM(means, cv2)
  #head(dm.stat)
  
  DM_HVG = names(dm.stat)[which(dm.stat>0.3)]
  #DM_HVG = DM_HVG[which()]
  #dev.off()
  
  sce.HVG.Brenneck = sce[rownames(sce)%in%Brennecke_HVG, ]
  sce.HVG.DM = sce[rownames(sce)%in%DM_HVG, ]
  
  save(sce, sce.HVG.Brenneck, sce.HVG.DM, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  #HVG_genes <- Brennecke_HVG$Gene
  #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  
  sce.sels = sce.HVG.Brenneck
  #sce.sels = sce.HVG.DM
  sce = sce.sels
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  # remove features with duplicated names
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  
  plotPCA(
    sce,
    run_args = list(exprs_values = "logcounts"),
    colour_by = "total_counts",
    size_by = "total_features_by_counts"
  )
  
  ### test SC3 clustering method
  sce = sc3_estimate_k(sce)
  metadata(sce)$sc3$k_estimation
  
  sce <- sc3(sce.sels, gene_filter = FALSE, ks = 2:30, biology = TRUE, n_cores = 6)
  #rowData(sce)$feature_symbol
  #sce <- sc3(sce, ks = 2, gene_filter = TRUE, biology = TRUE)
  #deng <- plotTSNE(deng, rand_seed = 1, return_SCE = TRUE)
  
  col_data <- colData(sce)
  head(col_data[ , grep("sc3_", colnames(col_data))])
  plotPCA(
    sce, 
    colour_by = "sc3_10_clusters", 
    size_by = "sc3_10_log2_outlier_score"
  )
  
  plotPCA(
    sce, 
    colour_by = "sc3_3_clusters",
    shape_by = "celltypes",
    size_by = "total_features_by_counts"
  )
  
  row_data <- rowData(sce)
  head(row_data[ , grep("sc3_", colnames(row_data))])
  
  plotFeatureData(
    sce, 
    aes(
      x = sc3_3_markers_clusts, 
      y = sc3_3_markers_auroc, 
      colour = sc3_3_markers_padj
    )
  )
  
  set.seed(1)
  plotTSNE(sce, colour_by="sc3_6_clusters", size_by = "total_features_by_counts",
           run_args = list(perplexity=20)) 
  
  #sc3_plot_consensus(sce, k = 3)
  sc3_plot_consensus(
    sce, k = 3, 
    show_pdata = c(
      "celltypes", 
      "sc3_3_clusters", 
      "log10_total_features_by_counts",
      "log10_total_counts",
      
      "sc3_3_log2_outlier_score"
    )
  )
  
  sc3_plot_cluster_stability(sce, k = 6)
  
  sc3_plot_de_genes(sce, k = 3, p.val = 0.2)
  
  sc3_plot_markers(sce, k = 3)
  
}




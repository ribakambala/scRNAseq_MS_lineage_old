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
# Convert Ensem ID to gene names
# a) keep only rows for mapped transripts
# b) map the gene names to the gene symbols
##########################################
convertGeneNames.forCountTable = function(aa, discard.gene.with.zero.reads = TRUE)
{
  
  aa = aa[ grep("^__", aa$gene, invert = TRUE), ] ## remove features that are not well aligned
  
  ## load annotation and change the gene names
  annot = read.csv(file = "../../../../annotations/BioMart_WBcel235_noFilters.csv", header = TRUE)
  gg.Mt = annot$Gene.name[which(annot$Chromosome.scaffold.name=="MtDNA")]
  gg.ribo = annot$Gene.name[which(annot$Gene.type=="rRNA")]
  
  ggs = as.character(aa$gene);
  mm = match(aa$gene, annot$Gene.stable.ID)
  jj = which(!is.na(mm)==TRUE & annot$Gene.name[mm] != "");
  ggs[jj] = as.character(annot$Gene.name[mm[jj]]);
  
  counts = aa
  rownames(counts) <- aa[, 1]
  counts <- as.matrix(counts[,-1])
  
  if(discard.gene.with.zero.reads){
    cat("-- discard genes with zero reads --")
    ss = apply(counts, 1, sum)
    keep.genes = which(ss>0)
    counts = counts[keep.genes, ]
    ggs = ggs[keep.genes]
    ggs.unique = unique(ggs)
    
    rownames(counts) = ggs
  }
  
  return(counts)
}

find.particular.geneSet = function(geneSet = "Mt")
{
  annot = read.csv(file = "../../../../annotations/BioMart_WBcel235_noFilters.csv", header = TRUE)
  if(geneSet == "Mt"){
    return(annot$Gene.name[which(annot$Chromosome.scaffold.name=="MtDNA")])
  }else{
    if(geneSet == "ribo"){
      return(annot$Gene.name[which(annot$Gene.type=="rRNA")])
    }else{
      stop("no geneSet found for :", geneSet, "\n")
    }
  }
    
}

##########################################
# function for collapse technical replicates in the lane level
# 
##########################################
compare.techinical.replicates = function(design.tr, counts.tr, filter.cells = FALSE, filter.genes = TRUE, check.correlations = TRUE)
{
  sinfos.uniq = unique(design.tr$seqInfos)
  
  ## add some new features for design for quality controls
  design.tr$log10_Total = log10(design.tr$total_reads)
  design.tr$percent_rRNA = design.tr$rRNA / design.tr$Total
  
  ## make SCE object and remove genes with zero reads detected
  sce <- SingleCellExperiment(assays = list(counts = counts.tr), 
                              colData = as.data.frame(design.tr), 
                              rowData = data.frame(gene_names = rownames(counts.tr), feature_symbol = rownames(counts.tr)))
  
  #is.spike <- grepl("^ERCC", rownames(sce))
  is.mito <- rownames(sce) %in% gg.Mt;
  is.ribo <- rownames(sce) %in% gg.ribo;
  summary(is.mito)
  summary(is.ribo)
  
  sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito, Ribo=is.ribo))
  
  # check some general statistics for each request and lane
  #cols = c(rep('gray', 2), 'red', 'darkblue', 'darkred', 'blue', 'red')
  ps1 = plotColData(sce, y = "log10_Total", x = "seqInfos") + ggtitle("total nb of reads")
  ps2 = plotColData(sce, y="uniquely_mapped_percent", x="seqInfos") + ggtitle("% of uniquely mapped ")
  ps3 = plotColData(sce, y="percent_assigned", x="seqInfos") + ggtitle("% of assigned")
  ps4 = plotColData(sce, y="pct_counts_Ribo", x="seqInfos") + ggtitle("% of rRNA contamination")
  ps5 = plotColData(sce, y="pct_counts_Mt", x="seqInfos") + ggtitle("% of Mt")
  
  ps6 = plotColData(sce, y="log10_total_counts", x="seqInfos") + ggtitle("total nb of reads mapped to transcripts")
  ps7 = plotColData(sce, y="total_features_by_counts", x="seqInfos") + ggtitle("total nb of genes")
  
  ps8 = plotColData(sce,
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
  if(filter.cells){
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
  }
  if(filter.genes){
    num.cells <- nexprs(sce, byrow=TRUE)
    ave.counts <- calcAverage(sce)
    genes.to.keep <- num.cells > 5 & ave.counts >= 10  & ave.counts <10^6  # detected in >= 2 cells, ave.counts >=5 but not too high
    summary(genes.to.keep)
    # remove mt and ribo genes
    genes.to.keep = genes.to.keep & ! rownames(sce) %in% gg.Mt & ! rownames(sce) %in% gg.ribo
    summary(genes.to.keep)
    
    sce <- sce[genes.to.keep, ]
    
  }
  
  
  for(kk in 1:8)
  {
    eval(parse(text = paste0("plot(ps", kk, ")")))
  }
  
  # select cells having technical replicates and normalize them  
  sce.qc = sce
  reducedDim(sce.qc) <- NULL
  endog_genes <- !rowData(sce.qc)$is_feature_control
  
  
  Methods.Normalization = c("cpm", "DESeq2", "scran")
  for(method in Methods.Normalization)
  {
    
    set.seed(1234567)
    cat('testing normalization method -- ', method, "\n")
    
    main = paste0(method, " normalization");
        
    if(method == "cpm") { ### cpm
      assay(sce.qc, "logcounts") <- log2(calculateCPM(sce.qc, use_size_factors = FALSE) + 1)
    }
    
    if(method == "DESeq2"){
      source("scRNAseq_functions.R")
      sizeFactors(sce.qc) = calculate.sizeFactors.DESeq2(counts(sce.qc))
      sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
    }
    
    if(method == "scran"){
      ## scran normalization (not working here, because negative scaling factor found)
      qclust <- quickCluster(sce.qc, min.size = 50)
      sce.qc <- computeSumFactors(sce.qc, clusters = qclust)
      sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
    }
    
    ps.norms = scater::plotPCA(
      sce.qc[endog_genes, ],
      run_args = list(exprs_values = "logcounts"), 
      size_by = "total_counts",
      #size_by = "total_features_by_counts",
      colour_by = "seqInfos"
    ) + ggtitle(paste0("PCA --", method))
    plot(ps.norms)
    
    ## check the correction of the same cells from different technical replicates
    bcs = unique(sce.qc$barcodes)
    
    if(check.correlations){
      correlations = c()
      
      #par(mfrow=c(2,3)) 
      for(n in 1:length(bcs))
      {
        # n = 1
        kk = which(sce.qc$barcodes == bcs[n])
        xx = as.data.frame(logcounts(sce.qc[, kk[match(sinfos.uniq, sce.qc$seqInfos[kk])]]))
        ss = apply(xx, 1, sum)
        xx = xx[ss>0, ]
        colnames(xx) = paste0(sinfos.uniq, ".", method)
        
        if(n <= 10) pairs(xx, lower.panel=NULL, upper.panel=panel.fitting)
        
        if(length(sinfos.uniq) == 2) correlations = c(correlations, cor(xx[,1], xx[, 2]))
        if(length(sinfos.uniq) == 3) correlations = rbind(correlations, c(cor(xx[, 1], xx[, 2]), cor(xx[, 1], xx[, 3]), cor(xx[, 2], xx[, 3])))
      }
      
      if(length(sinfos.uniq) == 2) {
        names(correlations) = paste0("cor_", sinfos.uniq, collapse = "_")
        hist(correlations, breaks = 10)
      }
      if(length(sinfos.uniq) == 3) {
        colnames(correlations) = paste0("cor_", sinfos.uniq[c(1,2,3)], sinfos.uniq[c(2,3,1)])
        pairs(correlations, lower.panel=NULL, upper.panel=panel.fitting)
      }
      #colnames(correlations) = c('rep0.vs.hiseq.rep1', 'rep0.vs.hiseq.rep2', 'hiseq.rep1.vs.hiseq.rep_2')
    }
    
  }
  
}

merge.techinical.replicates = function(design, counts, 
                                       sampleInfos.techRep = list(c("R7130_HHG5KBGX9_1", "R7130_HLWTCBGX9_1"),
                                                                  c("R7130_HHGHNBGX9_1", "R7130_CCYTEANXX_4", "R7133_CD2GTANXX_5")))
{
  
  for(n in 1:length(sampleInfos.techRep))
  {
    # check if technical replicates can be found
    techRep = sampleInfos.techRep[[n]]
    mm = match(techRep, unique(design$seqInfos))
    if(any(is.na(mm))) stop("Missed technical replicates : ", sampleInfos.techRep[which(is.na(mm))], "\n")
    
    sels = match(design$seqInfos, techRep)
    sels = which(!is.na(sels))
    
    design.keep = design[-sels, ]
    counts.keep = counts[, -sels]
    
    design.sels = design[sels, ]
    counts.sels = counts[, sels]
    
    cat("-- start to compare technical replicates for lanes :", techRep, "\n")
    
    compare.techinical.replicates(design.tr=design.sels, counts.tr=counts.sels)
    
    cat("-- start to merge technical replicates for lanes :", techRep, "\n")
    
    bcs = unique(design.sels$barcodes)
    
    design.merged = design.sels[match(bcs, design.sels$barcodes), ]
    design.merged$seqInfos = paste0(techRep[1], "_merged")
    
    # double chekc the barcode order
    for(index.bc in 1:length(bcs)){
      if(bcs[index.bc] != design.merged$barcodes[index.bc]) {
        cat("barcode order is wrong \n")
      }
    }
    
    merged = matrix(0, nrow = nrow(counts.sels), ncol = length(bcs))
    
    for(m in 1:length(techRep))
    {
      jj.trep = which(design.sels$seqInfos == techRep[m])
      counts.trep = counts.sels[, jj.trep]
      colnames(counts.trep) = design.sels$barcodes[jj.trep]
      counts.trep = counts.trep[, match(bcs, colnames(counts.trep))]
      counts.trep = as.matrix(counts.trep)
      counts.trep[which(is.na(counts.trep))] = 0
      merged = merged + counts.trep
    }
    colnames(merged) = design.merged$samples
    
    design = rbind(design.keep, design.merged)
    counts = cbind(counts.keep, merged)
  }
  
  return(list(design = design, counts = counts))
  
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

panel.fitting = function (x, y, bg = NA, pch = par("pch"), cex = 0.3, col='black') 
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

test.normalization = function(sce, Methods.Normalization = c("cpm", "DESeq2", "scran"), min.size = 100)
{
  #Methods.Normalization = "DESeq2" 
  
  for(method in Methods.Normalization)
  {
    sce.qc = sce
    set.seed(1234567)
    
    cat('test normalization method -- ', method, "\n")
    main = paste0(method);
    
    if(method == "raw") { # raw log counts
      assay(sce.qc, "logcounts") <- log2(counts(sce.qc) + 1)
    }
    
    if(method == "cpm") { ### cpm
      assay(sce.qc, "logcounts") <- log2(calculateCPM(sce.qc, use_size_factors = FALSE) + 1)
    }
    if(method == "UQ"){
      logcounts(sce.qc) <- log2(cal_uq_Hemberg(counts(sce.qc)) + 1)
    }
    
    if(method == "DESeq2"){
      sizeFactors(sce.qc) = calculate.sizeFactors.DESeq2(counts(sce.qc))
      sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
    }
    if(method == "downsample") {
      assay(sce.qc, "logcounts") <- log2(Down_Sample_Matrix(counts(sce.qc)) + 1)
    }
    
    if(method == "scran"){
      ## scran normalization (not working here, because negative scaling factor found)
      qclust <- quickCluster(sce.qc, min.size = min.size,  method = 'igraph')
      sce.qc <- computeSumFactors(sce.qc, clusters = qclust)
      sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
    }
    
    if(method == "TMM"|method == "DESeq2"|method == "UQ"|method == "scran"){
      summary(sizeFactors(sce.qc))
      range(sizeFactors(sce.qc))
      
      plot(sce.qc$total_counts/1e6, sizeFactors(sce.qc), log="xy", main = paste0(method), 
           xlab="Library size (millions)", ylab="Size factor",
           pch=16)
      #legend("bottomright", col=c("black"), pch=16, cex=1.2, legend = "size factor from scran vs total library size")
    }
    
    p1 = scater::plotPCA(
      sce.qc[endog_genes, ],
      run_args = list(exprs_values = "logcounts"), 
      size_by = "total_counts",
      #size_by = "total_features_by_counts",
      colour_by = "seqInfos"
    ) + ggtitle(paste0("PCA -- ", main))
    
    param.perplexity = 10;
    p2 = plotTSNE(
      sce.qc[endog_genes, ],
      run_args = list(exprs_values = "logcounts", perplexity = param.perplexity), 
      size_by = "total_counts",
      #size_by = "total_features_by_counts",
      colour_by = "seqInfos"  
    ) + ggtitle(paste0("tSNE - perplexity = ", param.perplexity, "--", main))
    
    p3 = plotUMAP(
      sce.qc[endog_genes, ],
      run_args = list(exprs_values = "logcounts"), 
      size_by = "total_counts",
      #size_by = "total_features_by_counts",
      colour_by = "seqInfos"
    ) + ggtitle(paste0("UMAP -- ", main))
    
    plot(p1); plot(p2); plot(p3)
  }
  
}

########################################################
########################################################
# Section : funcitons for cell cycle scoring and correction
#  
########################################################
########################################################
find.cc.markers.homologues = function()
{
  detach("package:Seurat", unload=TRUE)
  require(Seurat)
  #s.genes = c("cdk-4", "evl-18") # from GO:1901987 http://amigo1.geneontology.org/cgi-bin/amigo/term-assoc.cgi?term=GO:1902808&speciesdb=all&taxid=6239
  #g2m.genes = xx # GO:1902751
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  homologues = read.delim("../../../../annotations/cellCycle_genes_worm/BioMart_worm_human_homologe.txt", sep = "\t",
                          header = TRUE)
  #homologues = homologues[which(homologues$Human.orthology.confidence..0.low..1.high.==1), ]
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  s.genes = homologues$Gene.name[match(s.genes, homologues$Human.gene.name)]
  g2m.genes = homologues$Gene.name[match(g2m.genes, homologues$Human.gene.name)]
  s.genes = s.genes[which(!is.na(s.genes))]
  g2m.genes = g2m.genes[which(!is.na(g2m.genes))]
  
  c3.genes = list(s.genes = s.genes, g2m.genes = g2m.genes)
  return(c3.genes)
}

find.cc.markers.GO = function()
{
  s.genes = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_dnasynthesis.txt", 
                       sep = "\t", header = FALSE)
  gm1 = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_G2Mphase.txt", 
                   sep = "\t", header = FALSE)
  gm2 = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_Mphase.txt", 
                   sep = "\t", header = FALSE)
  gm3 = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_G2phase.txt", 
                   sep = "\t", header = FALSE)
  s.genes = unique(s.genes[, 3])
  g2m.genes = c(unique(as.character(gm1[,3])), unique(as.character(gm2[, 3])), unique(as.character(gm3[, 3])))
  
  c3.genes = list(s.genes = s.genes, g2m.genes = g2m.genes)
  return(c3.genes)
}
  
find.cellcycle.markers = function(list.sel = 'homologues')
{
  ##########################################
  # different strategies to find cell cycle markers of c. elegans for Seurat
  # 1st method): using homologue between human and c. elegans
  # 2nd method): find a list of 
  ##########################################
  if(list.sel == 'homologues') c3.genes = find.cc.markers.homologues() 
  if(list.sel == "curated") c3.genes = find.cc.markers.GO()
  s.genes = c3.genes$s.genes
  g2m.genes = c3.genes$g2m.genes
  
  # manually add genes from wormbook http://www.wormbook.org/chapters/www_cellcyclereguln/cellcyclereguln.html
  #s.genes = c(as.character(s.genes), c("cye-1")
  s.genes = unique(c(as.character(s.genes), c('cye-1', 'cya-1', 'evl-18')))
  g2m.genes = unique(c(as.character(g2m.genes), c('cdk-1', 'mat-1', 'mat-2', 'mat-3', 'emb-27', 'emb-30', 'mdf-1', 'san-1')))
  
  c3.genes = list(s.genes = s.genes, g2m.genes = g2m.genes)
  return(c3.genes)
  
  # test the code from https://github.com/hbc/macrae_ghazizadeh_zebrafish_heart_sc/blob/master/seurat_cluster/seurat_cluster_adapted_WT.Rmd
  # unfornately it does not work, because several pacakges can not be properly installed
  Test.query.cellCycle.markers = FALSE
  if(Test.query.cellCycle.markers){
    require(plotly)
    require(remotes)
    annot <- basejump::annotable("Danio rerio") %>% 
      dplyr::select(c(ensgene,symbol)) %>% 
      dplyr::mutate(symbol = toupper(symbol)) 
    cell_cycle_markers <- bcbioSingleCell::cellCycleMarkers[[camel("mus musculus")]] %>% 
      dplyr::mutate(symbol = toupper(symbol)) %>% dplyr::inner_join(annot,by = "symbol") %>% 
      dplyr::select(-c(ensgene.x)) %>% dplyr::rename(ensgene = ensgene.y)
    stopifnot(is.data.frame(cell_cycle_markers))
    markdownHeader("S phase markers", level = 3)
    s_genes <- cell_cycle_markers %>%
      filter(phase == "S") %>%
      pull("ensgene")
    print(s_genes)
    markdownHeader("G2/M phase markers", level = 3)
    g2m_genes <- cell_cycle_markers %>%
      filter(phase == "G2/M") %>%
      pull("ensgene")
    print(g2m_genes)
    saveData(cell_cycle_markers, s_genes, g2m_genes, dir = data_dir)
  }
  
}

cellCycle.correction = function(sce, method = "seurat")
{
  
  if(method == 'scran'){
    set.seed(100)
    library(scran)
    mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                    package="scran"))
    assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
    plot(assignments$score$G1, assignments$score$G2M, 
         xlab="G1 score", ylab="G2/M score", pch=16)
    sce$phases <- assignments$phases
    table(sce$phases)
  }
  
  if(method == "seurat"){
    library(scater)
    # install loomR from GitHub using the remotes package remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
    library(loomR)
    library(Seurat)
    seurat = as.Seurat(sce, counts = "counts", data = "logcounts")
    detach("package:scater", unload=TRUE)
    
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 1000)
    
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(seurat), 25)
    plot1 <- VariableFeaturePlot(seurat)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    CombinePlots(plots = list(plot1, plot2))
    
    seurat <- ScaleData(seurat, features = rownames(seurat), model.use = "linear")
    seurat <- RunPCA(seurat, features = VariableFeatures(seurat), ndims.print = 6:10, nfeatures.print = 10)
    
    DimPlot(seurat, reduction = "pca")
    DimHeatmap(seurat, dims = c(8, 10))
    
    source("scRNAseq_functions.R")
    c3.genes = find.cellcycle.markers(list.sel = "homologues")
    
    s.genes <- c3.genes$s.genes
    g2m.genes <- c3.genes$g2m.genes 
    s.genes = s.genes[which(!is.na(match(s.genes, rownames(seurat))))]
    g2m.genes = g2m.genes[which(!is.na(match(g2m.genes, rownames(seurat))))]
    
    seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    
    # view cell cycle scores and phase assignments
    # head(seurat[[]])
    
    RidgePlot(seurat, features = c("cdk-1", "cdk-4", "cyd-1", "cye-1", "cya-1", "wee-1.3"), ncol = 2)
    seurat <- RunPCA(seurat, features = VariableFeatures(seurat))
    p0 = DimPlot(seurat)
    
    seurat <- RunPCA(seurat, features = c(as.character(s.genes), as.character(g2m.genes)))
    DimPlot(seurat)
    
    # regress out the cell cycle
    seurat1 <- ScaleData(seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat))
        
    seurat1 <- RunPCA(seurat1, features = VariableFeatures(seurat1), nfeatures.print = 10)
    p1 = DimPlot(seurat1)
    
    # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
    seurat1 <- RunPCA(seurat1, features = c(s.genes, g2m.genes))
    DimPlot(seurat1)
    
    # regressing out the difference between the G2M and S phase scores
    seurat$CC.Difference <- seurat$S.Score - seurat$G2M.Score
    seurat2 <- ScaleData(seurat, vars.to.regress = "CC.Difference", features = rownames(seurat))
    
    # cell cycle effects strongly mitigated in PCA
    seurat2 <- RunPCA(seurat2, features = VariableFeatures(seurat2), nfeatures.print = 10)
    p2 = DimPlot(seurat2)
    
     
    # when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1
    # cells however, within actively proliferating cells, G2M and S phase cells group together
    seurat2 <- RunPCA(seurat2, features = c(s.genes, g2m.genes))
    DimPlot(seurat2)
    
    # save cell cycle scoring and corrected matrix
    library(scater)
    sce$S.Score = seurat$S.Score
    sce$G2M.Score = seurat$G2M.Score
    sce$Phase = seurat$Phase
    #sce$Phase.GO = seurat$old.ident
    sce$CC.Difference = seurat$CC.Difference
    
    xx = as.data.frame(seurat@assays$RNA@scale.data); rownames(xx) = rownames(sce)
    assay(sce, "logcounts_seurat") <- xx
    xx = as.data.frame(seurat1@assays$RNA@scale.data); rownames(xx) = rownames(sce)
    assay(sce, "logcounts_seurat_cellcycleCorrected") <- xx
    xx = as.data.frame(seurat2@assays$RNA@scale.data); rownames(xx) = rownames(sce)
    assay(sce, "logcounts_seurat_SG2MCorrected") <- xx
    
    save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected.Rdata'))
    
    return(sce)
    
    ##########################################
    # a peek of source code for cellCycle.scoring function from Seurat
    # https://github.com/ChristophH/in-lineage/blob/master/R/lib.R
    ##########################################
    Example.for.Seurat = FALSE
    if(Example.for.Seurat){
      library('Matrix')
      library('parallel')
      library('MASS')
      library('diffusionMap')
      library('FNN')
      library('igraph')
      library('princurve')
      library('ggplot2')
      library('inline')
      library('gplots')
      
      # for cell cycle score
      get.bg.lists <- function(goi, N, expr.bin) {
        res <- list()
        goi.bin.tab <- table(expr.bin[goi])
        for (i in 1:N) {
          res[[i]] <- unlist(lapply(names(goi.bin.tab), function(b) {
            sel <- which(expr.bin == as.numeric(b) & !(names(expr.bin) %in% goi))
            sample(names(expr.bin)[sel], goi.bin.tab[b])
          }))
        }
        return(res)
      }
      
      enr.score <- function(expr, goi, bg.lst) {
        goi.mean <- apply(expr[goi, ], 2, mean)
        bg.mean <- sapply(1:length(bg.lst), function(i) apply(expr[bg.lst[[i]], ], 2, mean))
        return((goi.mean - apply(bg.mean, 1, mean)) / apply(bg.mean, 1, sd))
      }
      
      get.cc.score <- function(cm, N=100, seed=42, 
                               s.gene.file='./annotation/s_genes.txt',
                               g2m.gene.file='./annotation/g2m_genes.txt')
      {
        set.seed(seed)
        cat('get.cc.score, ')
        cat('number of random background gene sets set to', N, '\n')
        
        min.cells <- 5
        
        cells.mols <- apply(cm, 2, sum)
        gene.cells <- apply(cm>0, 1, sum)
        cm <- cm[gene.cells >= min.cells, ]
        
        gene.mean <- apply(cm, 1, mean)
        
        breaks <- unique(quantile(log10(gene.mean), probs = seq(0,1, length.out = 50)))
        gene.bin <- cut(log10(gene.mean), breaks = breaks, labels = FALSE)
        names(gene.bin) <- rownames(cm)
        gene.bin[is.na(gene.bin)] <- 0
        
        regev.s.genes <- read.table(file=s.gene.file, header=FALSE, stringsAsFactors=FALSE)$V1
        regev.g2m.genes <- read.table(file=g2m.gene.file, header=FALSE, stringsAsFactors=FALSE)$V1
        
        goi.lst <- list('S'=rownames(cm)[!is.na(match(toupper(rownames(cm)), regev.s.genes))],
                        'G2M'=rownames(cm)[!is.na(match(toupper(rownames(cm)), regev.g2m.genes))])
        
        n <- min(40, min(sapply(goi.lst, length)))
        goi.lst <- lapply(goi.lst, function(x) x[order(gene.mean[x], decreasing = TRUE)[1:n]])
        
        bg.lst <- list('S'=get.bg.lists(goi.lst[['S']], N, gene.bin),
                       'G2M'=get.bg.lists(goi.lst[['G2M']], N, gene.bin))
        
        all.genes <- sort(unique(c(unlist(goi.lst, use.names=FALSE), unlist(bg.lst, use.names=FALSE))))
        
        expr <- log10(cm[all.genes, ]+1)
        
        s.score <- enr.score(expr, goi.lst[['S']], bg.lst[['S']])
        g2m.score <- enr.score(expr, goi.lst[['G2M']], bg.lst[['G2M']])
        
        phase <- as.numeric(g2m.score > 2 & s.score <= 2)
        phase[g2m.score <= 2 & s.score > 2] <- -1
        
        return(data.frame(score=s.score-g2m.score, s.score, g2m.score, phase))
      }
    }
    
  }
  
  if(method == "scLVM"){
    ##########################################
    # select the python verson to use for Rstudio
    # https://cran.r-project.org/web/packages/reticulate/vignettes/versions.html
    # still does not work at the end
    # we change the strategy: prepare the tables and run scLVM in the conda version
    ##########################################
    # system("python --version")
    # system("which python")
    # 
    # library(reticulate)
    # use_python("/Users/jiwang/anaconda3/envs/scLVM/bin/python")
    # #use_condaenv(condaenv = "scLVM", conda = "/Users/jiwang/anaconda3/condabin/conda", required = TRUE)
    # py_config()
    # 
    # system("python --version")
    # system("which python")
    # 
    # Sys.setenv(PATH = paste("/Users/jiwang/anaconda3/envs/scLVM/bin", Sys.getenv("PATH"),sep=":"))
    
    #install.packages("rPython", type = "source")
    #install.packages("/Users/jiwang/src_packages/scLVM_0.99.3.tar.gz", repos = NULL, type="source")
    
    ##########################################
    # example code from scLVM R tutorial
    # https://github.com/PMBio/scLVM/blob/master/R/tutorials/scLVM_vignette.Rmd
    ##########################################
    library(rPython)
    library(genefilter)
    library(statmod)
    require(ggplot2)
    library(gplots)
    require(DESeq2)
    library(scLVM)
    
    #limix_path = '/Users/jiwang/anaconda2/envs/scLVM/bin/python'
    #configLimix(limix_path)
    
    data(data_Tcells)
    help(data_Tcells)

    #dataMouse[ 1:5, 1:4 ]
    geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
      substr( rownames(dataMouse), 1, 4 ) ] )
    #2. calculate normalisation for counts
    countsMmus <- dataMouse[ which( geneTypes=="ENSM" ), ]
    countsERCC <- dataMouse[ which( geneTypes=="ERCC" ), ]
    lengthsMmus <- dataMouse[ which( geneTypes=="ENSM" ), 1 ]
    lengthsERCC <- dataMouse[ which( geneTypes=="ERCC" ), 1 ]
    sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
    sfMmus <- sfERCC #also use ERCC size factor for endogenous genes
    #normalise read counts
    nCountsERCC <- t( t(countsERCC) / sfERCC )
    nCountsMmus <- t( t(countsMmus) / sfERCC )
    
    countsMmus = counts(sce)
    sfERCC = estimateSizeFactorsForMatrix(countsMmus)
    sfMmus <- sfERCC
    nCountsMmus = t( t(countsMmus) / sfERCC )
    #use spike in to find tehcnical noise. 
    # If no spike-ins are available, we can also use the endogenous read counts for fitting the mean-CV2 relation using a log-linear fit in the log-space.
    # Alternatively, we can fit the mean-variance relationship in the log-space using local 2nd order polynomial regression (loess).
    #techNoise = fitTechnicalNoise(nCountsMmus,nCountsERCC=nCountsERCC, fit_type = 'counts')  
    techNoiseLogFit = fitTechnicalNoise(nCountsMmus, fit_type = 'log', use_ERCC = FALSE, plot=TRUE) 
    #techNoiseLogVarFit = fitTechnicalNoise(nCountsMmus, fit_type = 'logvar', use_ERCC = FALSE, plot=TRUE) 
    
    #call variable genes
    #is_het = getVariableGenes(nCountsMmus, techNoiseLogFit$fit, method = "fit", 
    #                          threshold = 0.1, fit_type="log",sfEndo=sfMmus, sfERCC=sfERCC)
    #table(is_het)
    
    #we an also do this for the other fits
    is_hetLog = getVariableGenes(nCountsMmus, techNoiseLogFit$fit, plot=TRUE)
    table(is_hetLog)
    #is_hetLogVar = getVariableGenes(nCountsMmus, techNoiseLogVarFit$fit, plot=TRUE)
    #table(is_hetLogVar)
    
    #get cell cycle genes from GO
    cc.genes = read.delim("../../../../annotations/cellCycle_genes_worm/GO_0007049_genes.txt", 
                         sep = "\t", header = FALSE)
    cc.genes = unique(cc.genes[,3])
    
    #rename a few variables
    Y = t(log10(nCountsMmus+1)) #normalised trandformed read counts
    genes_het_bool = as.vector(is_hetLog) #variable genes
    #genes_het_bool[]
    geneID = rownames(nCountsMmus) #gene IDs
    tech_noise = as.vector(techNoiseLogFit$techNoiseLog) #technical noise
    ens_ids_cc <- cc.genes
    
    index.cc = match(cc.genes, geneID)
    index.cc = index.cc[which(!is.na(index.cc))]
    ##########################################
    # can not proceed anymore and save tables for python in conda
    ##########################################
    #construct and initialize new scLVM object
    sclvm = new("scLVM")
    sclvm = init(sclvm,Y=Y,tech_noise = tech_noise)
    
    # CellCycleARD = fitFactor(sclvm,geneSet = ens_ids_cc, k=20,use_ard = TRUE)
    
    write.table(Y, file = paste0(tabDir, "gene_expression_matrx_4scLVM.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(tech_noise, file = paste0(tabDir, "tech_noise_4scLVM.txt"), sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE )
    write.table(geneID, file =paste0(tabDir, "geneNames_4scLVM.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE )
    write.table(which(genes_het_bool==TRUE), file =paste0(tabDir, "index_hetgenes_4scLVM.txt"), sep = "\t", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    write.table(index.cc, file =paste0(tabDir, "index_ccgenes_4scLVM.txt"), sep = "\t", 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  if(method == "ccRemover"){
    # see examplel from https://cran.r-project.org/web/packages/ccRemover/vignettes/ccRemover_tutorial.html
    # this method is pretty slow and probably risky (to test), because simply PCA was done and 
    # then compare cell cycle genes were compares with control genes and significantly different PCs were removed
    require(ccRemover)
    t.cell_data = logcounts(sce)
    head(t.cell_data[,1:5])
    
    summary(apply(t.cell_data,1, mean))
    mean_gene_exp <- rowMeans(t.cell_data)
    t_cell_data_cen <- t.cell_data - mean_gene_exp
    summary(apply(t_cell_data_cen,1,mean))
    
    gene_names <- rownames(t_cell_data_cen)
    # cell_cycle_gene_indices <- gene_indexer(gene_names, species = "mouse", 
    #                                         name_type = "symbols" )
    # length(cell_cycle_gene_indices)
    # if_cc <- rep(FALSE,nrow(t_cell_data_cen)) 
    # if_cc[cell_cycle_gene_indices] <- TRUE
    # summary(if_cc)
    
    cc.genes = read.delim("../../../../annotations/cellCycle_genes_worm/GO_0007049_genes.txt", header = FALSE)
    cc.genes = unique(cc.genes$V3)
    cc.index = match(cc.genes, gene_names)
    cc.index = cc.index[which(!is.na(cc.index))]

    if_cc <- rep(FALSE,nrow(t_cell_data_cen))
    if_cc[cc.index] <- TRUE
    summary(if_cc)
    
    dat <- list(x=t_cell_data_cen, if_cc=if_cc)
    xhat <- ccRemover(dat, bar=TRUE, max_it = 6, nboot = 100)
    
    xhat <- xhat + mean_gene_exp
    
    pca1 = prcomp(t(t.cell_data[if_cc,]), scale. = TRUE)
    pca2 = prcomp(t(xhat[if_cc,]), scale. = TRUE)
    par(mfrow = c(1, 2))
    plot(pca1$x[, c(2:3)])
    plot(pca2$x[, c(1:2)])
    
  }
  
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




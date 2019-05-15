##################################################
##################################################
## Project: transcriptional priming mechanism by tbx in C.elegans
## Script purpose: analysis the single-cell RNA-seq data
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Feb 19 14:43:38 2018
##################################################
##################################################
version.DATA = 'R5399_scRNA_v1'
version.analysis =  paste0(version.DATA, '_2018_02_19')

design.file = "../exp_design/R5399_NGS_Analysis_Design.csv"
dataDir = paste0("../data/", version.DATA)

resDir = paste0("../results/", version.analysis)
tabDir = paste0("../results/", version.analysis, "/tables")
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}

##################################################
##################################################
## Section: import sample design and data and prepare table 
##################################################
##################################################
design = read.csv(design.file, header = TRUE, as.is = c(1))
colnames(design) = c("sampleID", "sampleInfo", "celltypes")
design = design[which(!is.na(design$sampleID)), ]
design$celltypes = sapply(design$celltypes, function(x) gsub("\\(", "", as.character(x)))
design$celltypes = sapply(design$celltypes, function(x) gsub("\\)", "", as.character(x)))
design$celltypes = sapply(design$celltypes, function(x) gsub(" ", ".", as.character(x)))

#conditions = sapply(design$sampleInfo, function(x) paste0(rev(rev(unlist(strsplit(as.character(x), "_")))[-c(1:3)]), collapse="_"))
#times =  sapply(design$sampleInfo, function(x) rev(unlist(strsplit(as.character(x), "_")))[3])

design = data.frame(design[, c(1, 3, 2)],stringsAsFactors = FALSE)
#design$conditions = sapply(design$conditions, function(x) gsub("_", "-", x))

xlist <-list.files(path=dataDir, pattern = "*.txt", full.names = TRUE)
#ylist = list.files(path=path, pattern = "Information_design*", full.names = TRUE)
aa = NULL
for(n in 1:length(xlist))
{
  cat(n, '\t')
  cat(xlist[n], '\n')
  
  if(n==1){
    aa = read.table(xlist[n], header = FALSE);
    colnames(aa) = c('RefseqID', basename(xlist[n]))
    aa = data.frame(aa, stringsAsFactors = FALSE)
  }else{
    test = read.table(xlist[n], header = FALSE);
    colnames(test) = c('RefseqID',  basename(xlist[n]))
    mm = match(aa$RefseqID, test$RefseqID)
    aa = data.frame(aa, test[mm, -grep('RefseqID', colnames(test))])
  }
}
colnames(aa)[1] = 'gene';
colnames(aa)[-1] = sapply(xlist, basename)

index = c()
for(n in 1:nrow(design))
{
  jj = grep(design$sampleID[n], colnames(aa))
  if(length(jj)==1){index = c(index, jj)
  }else{cat("NOT FOUND sample for ", design$sampleID[n], "\n")}
}

aa = data.frame(aa$gene, aa[, index], stringsAsFactors = FALSE)
colnames(aa)[-1] = apply(design[, c(2,1)], 1, paste0, collapse="_")
colnames(aa)[1] = 'gene';
save(design, aa, file=paste0('../Rdata/', version.DATA, '_RAW_Read_Counts_RNA_seq.Rdata'))

##################################################
##################################################
## Section: Quality control and clean the count table for scRNA-seq
##################################################
##################################################
load(file=paste0('../Rdata/', version.DATA, '_RAW_Read_Counts_RNA_seq.Rdata'))
aa = aa[ grep("^__", aa$gene, invert = TRUE), ] ## remove features that are not well aligned

## load annotation and change the gene names
load(file = "../../../../annotations/BioMart_WBcel235.Rdata")
gg.Mt = annot$Gene.name[which(annot$Chromosome.scaffold.name=="MtDNA")]

ggs = as.character(aa$gene);
mm = match(aa$gene, annot$Gene.stable.ID)
jj = which(!is.na(mm)==TRUE & annot$Gene.name[mm] != "");
ggs[jj] = as.character(annot$Gene.name[mm[jj]]);

counts = aa
rownames(counts) <- aa[, 1]
counts <- as.matrix(counts[,-1])
#celltype <- design$celltypes

ss = apply(counts, 1, sum)
keep.genes = which(ss>0)
counts = counts[keep.genes, ]
ggs = ggs[keep.genes]
ggs.unique = unique(ggs)

rownames(counts) = ggs

## import packages for the QC and table cleaning
library(SingleCellExperiment)
library(scater)
require(scran)
library(limma)
options(stringsAsFactors = FALSE)

## make SCE object and remove genes with zero reads detected
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = design, 
                            rowData = data.frame(gene_names = rownames(counts), feature_symbol = rownames(counts)))

write.csv(counts(sce), file=paste0(tabDir, "/scRNAseq_raw_readCounts.csv"), row.names=TRUE)

keep_feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_feature, ]


#is.spike <- grepl("^ERCC", rownames(sce))
is.mito <- rownames(sce) %in% gg.Mt;
summary(is.mito)

sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito))
head(colnames(colData(sce)))

####################
## filter cells with low quality 
####################
pdfname = paste0(resDir, "/scRNAseq_QCs_cells_filterting_", version.analysis, ".pdf")
pdf(pdfname, width=10, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
#par(mfrow=c(1, 1))
# par(mfcol=c(1, 1))

par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$log10_total_counts, xlab="Library sizes (in log10)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
abline(v=c(6, 5, 4), col="red", lwd=2.0)

hist(sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
#hist(sce$pct_counts_ERCC, xlab="ERCC proportion (%)", 
#     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

par(mfrow=c(1, 1))
plot(sce$total_counts, sce$total_features, xlab="Library size", ylab='Nb of expressed genes', log='xy')
abline(v=c(10^4, 10^5, 10^6), col='red', lwd=2.0)
abline(h=c(100, 200, 1000, 2000), col='red', lwd=2.0)

plotPhenoData(sce, 
              aes_string(
                x = "log10_total_counts",
                y = "log10_total_features",
                colour = "celltypes",
                size = 12
              )
)

plotPhenoData(sce, 
              aes_string(
    x = "total_features",
    y = "pct_counts_Mt",
    colour = "celltypes",
    size = 12
  )
)
## here we are using the 50,000 for library size and 100 expressed genes as thresholds
threshod.total.counts.per.cell = 10^5
threshod.nb.detected.genes.per.cell = 100;
#libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
#feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
#libsize.drop = sce$total_counts<5*10^4
#feature.drop = sce$total_features < 200
filter_by_total_counts <- (sce$total_counts > threshod.total.counts.per.cell)
table(filter_by_total_counts)
filter_by_expr_features <- (sce$total_features > threshod.nb.detected.genes.per.cell)
table(filter_by_expr_features)
filter_by_MT = sce$pct_counts_Mt < 5
table(filter_by_MT)

sce$use <- ( 
    filter_by_expr_features & # sufficient features (genes)
    filter_by_total_counts & # sufficient molecules counted
    # filter_by_ERCC & # sufficient endogenous RNA
    filter_by_MT # remove cells with unusual number of reads in MT genes
)
table(sce$use)

# check the cell filtering in PCA plot
sce <- plotPCA(
  sce,
  size_by = "total_features", 
  shape_by = "use",
  pca_data_input = "pdata",
  detect_outliers = TRUE,
  return_SCE = TRUE
)
table(sce$outlier)

#fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
#plotPCA(sce, pca_data_input="pdata") + fontsize

## comparison between manual filtering and automatic ouliter filtering
Manual.vs.outlier.filtering = TRUE
if(Manual.vs.outlier.filtering)
{
  auto <- colnames(sce)[sce$outlier]
  man <- colnames(sce)[!sce$use]
  venn.diag <- vennCounts(
    cbind(colnames(sce) %in% auto,
          colnames(sce) %in% man)
  )
  vennDiagram(
    venn.diag,
    names = c("Automatic", "Manual"),
    circle.col = c("blue", "green")
  )
}


dev.off()

sce = sce[, sce$use]
save(sce, file=paste0('../Rdata/', version.DATA, '_QCed_cells_filtered_SCE.Rdata'))

####################
## filter lowly expressed (and probably too highly expressed genes)
####################
pdfname = paste0(resDir, "/scRNAseq_QCs_genes_filterting_", version.analysis, ".pdf")
pdf(pdfname, width=10, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

#plotQC(reads, type = "highest-expression", n=20)
fontsize <- theme(axis.text=element_text(size=16), axis.title=element_text(size=16))
plotQC(sce, type = "highest-expression", n=20) + fontsize

ave.counts <- calcAverage(sce)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
     xlab=expression(Log[10]~"average count"))

num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", 
              xlab=expression(Log[10]~"average count"))
dev.off()

genes.to.keep <- num.cells > 2 & ave.counts >= 1  & ave.counts <10^5  # detected in >= 2 cells, ave.counts >=5 but not too high
summary(genes.to.keep)
sce <- sce[genes.to.keep, ]

save(sce, file=paste0('../Rdata/', version.DATA, '_QCed_cells_genes_filtered_SCE.Rdata'))

##################################
##### data visulization and identify confounding factors
##################################
Data.visulization.exploration = FALSE
if(Data.visulization.exploration)
{
  load(file=paste0('../Rdata/', version.DATA, '_QCed_cells_genes_filtered_SCE.Rdata'))
  sce.qc = sce
  
  assay(sce, "logcounts_raw") <- log2(counts(sce) + 1)
  reducedDim(sce) <- NULL
  #reducedDim(sce.qc) <- NULL
  endog_genes <- !rowData(sce)$is_feature_control
  
  pdfname = paste0(resDir, "/scRNAseq_filtered_identifying_confounding_factors_", version.analysis, ".pdf")
  pdf(pdfname, width=10, height = 6)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  plotPCA(
    sce[endog_genes, ],
    exprs_values = "counts",
    colour_by = "celltypes",
    size_by = "total_features"
    #shape_by = "individual"
  ) + ggtitle("PCA plot using counts")
  
  plotPCA(
    sce[endog_genes, ],
    exprs_values = "logcounts_raw",
    colour_by = "celltypes",
    size_by = "total_features"
    #shape_by = "individual"
  ) + ggtitle("PCA plot using logcounts_raw")
  
  param.perplexity = 5;
  plotTSNE(
    sce[endog_genes, ],
    exprs_values = "logcounts_raw",
    perplexity = param.perplexity,
    colour_by = "celltypes",
    size_by = "total_features",
    rand_seed = 123456
  ) + ggtitle("t-SNE using logcounts_raw")
  
  ##### identify confounding factors
  plotQC(
    sce[endog_genes, ],
    type = "find-pcs",
    exprs_values = "logcounts_raw",
    variable = "total_features"
  )
  
  plotQC(
    sce[endog_genes, ],
    type = "expl",
    exprs_values = "logcounts_raw",
    variables = c(
      "total_features",
      "total_counts",
      "pct_counts_Mt"
    )
  )
  
  dev.off()

}

##################################
##### normalization
##################################
load(file=paste0('../Rdata/', version.DATA, '_QCed_cells_genes_filtered_SCE.Rdata'))
assay(sce, "logcounts_raw") <- log2(counts(sce) + 1)
reducedDim(sce) <- NULL
sce.qc = sce
endog_genes <- !rowData(sce)$is_feature_control

library(scran)
options(stringsAsFactors = FALSE)
set.seed(1234567)

pdfname = paste0(resDir, "/scRNAseq_filtered_normalization_", version.analysis, ".pdf")
pdf(pdfname, width=10, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

### raw counts in log scale
plotPCA(
  sce.qc[endog_genes, ],
  exprs_values = "logcounts_raw",
  colour_by = "celltypes",
  size_by = "total_features"
) + ggtitle("PCA plot using logcounts_raw")

### cpm
assay(sce, "logCPM") <- log2(calculateCPM(sce, use.size.factors = FALSE) + 1)
plotPCA(
  sce[endog_genes, ],
  exprs_values = "logCPM",
  colour_by = "celltypes",
  size_by = "total_features"
) + ggtitle("PCA plot using logCPM")

plotRLE(sce[endog_genes, ], 
  exprs_mats = list(Raw = "logcounts_raw", CPM = "logCPM"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "celltypes"
)
## upperquartile normalization
sce <- normaliseExprs(
  sce,
  method = "upperquartile", 
  feature_set = endog_genes,
  p = 0.99,
  return_log = TRUE,
  return_norm_as_exprs = TRUE
) 

plotPCA(
  sce[endog_genes, ],
  colour_by = "celltypes",
  size_by = "total_features"
) + ggtitle("PCA plot using log(normalized by upperquartile)")

## size factor normalization (RLE)
sce <- normaliseExprs(sce,
                      method = "RLE", 
                      feature_set = endog_genes,
                      return_log = TRUE,
                      return_norm_as_exprs = TRUE
)

plotPCA(
  sce[endog_genes, ],
  colour_by = "celltypes",
  size_by = "total_features"
) + ggtitle("PCA plot using log(normalized by RLE )")

## scran normalization (not working here, because negative scaling factor found)
#qclust <- quickCluster(sce.qc, min.size = 2)
#reads.qc <- computeSumFactors(reads.qc, sizes = 5, clusters = qclust)
#umi.qc <- normalize(umi.qc)
sce <- computeSumFactors(sce, sizes = seq(2, 8, 2))
summary(sizeFactors(sce))

if(any(sizeFactors(sce)<0)){
  cat("ERROR........")
  cat("NEGATIVE size factors FOUND !!! \n" )
  cat(sizeFactors(sce), "\n")
}else{
  plot(sizeFactors(sce), sce$total_counts/1e6, log="xy",
       ylab="Library size (millions)", xlab="Size factor")
  sce <- normalize(sce)
  
  plotPCA(
    sce[endog_genes, ],
    exprs_values = "normcounts",
    colour_by = "celltypes",
    size_by = "total_features"
  )
}

dev.off()

save(sce, file=paste0('../Rdata/', version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))

write.csv(logcounts(sce), file=paste0(tabDir, "/scRNAseq_QC_filtered_normalized_readCounts.csv"), row.names=TRUE)

## check clustering after normalization
plotPCA(
  sce[endog_genes, ],
  colour_by = "celltypes",
  size_by = "total_features"
) + ggtitle("PCA plot using log(normalized counts)")

param.perplexity = 6;
plotTSNE(
  sce[endog_genes, ],
  exprs_values = "logcounts",
  perplexity = param.perplexity,
  colour_by = "celltypes",
  size_by = "total_features",
  rand_seed = 123456
) + ggtitle("PCA plot using log(normalized counts)")

dev.off()


##################################################
##################################################
## Section: feature selection and dimension reduction) and clustering
##################################################
##################################################
library(scater)
library(scran)
library(scRNA.seq.funcs)
library(matrixStats)
library(M3Drop)
library(RColorBrewer)
library(SingleCellExperiment)
set.seed(1)
load(file=paste0('../Rdata/', version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))

var.fit.nospike <- trendVar(sce, parametric=TRUE, use.spikes=FALSE, assay.type="logcounts", span=0.2)
## the following funciton needs spike-in
#var.fit.nospike = technicalCV2(sce, spike.type=NULL, assay.type="normcounts", min.bio.disp=0.25)
#var.fit.nospike = improvedCV2(sce, use.spikes= FALSE, assay.type="normcounts",logged=FALSE, normalized=TRUE)
var.out.nospike <- decomposeVar(sce, var.fit.nospike)

#plot(var.fit.nospike$mean, var.fit.nospike$var)

pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_HVG_featureSelect_clustering_", version.analysis, ".pdf")
pdf(pdfname, width=10, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
par(mfrow=c(1, 1))
# par(mfcol=c(1, 1))

plot(var.out.nospike$mean, var.out.nospike$total, pch=16, cex=0.6, 
     xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=TRUE)
#points(var.out.nospike$mean[cur.spike], var.out.nospike$total[cur.spike], col="red", pch=16)

chosen.genes <- order(var.out.nospike$bio, decreasing=TRUE)[1:50]
plotExpression(sce, features=rownames(var.out.nospike)[chosen.genes]) + fontsize

sce <- denoisePCA(sce, technical=var.fit.nospike$trend) 
dim(reducedDim(sce, "PCA")) 

plotReducedDim(sce, use_dimred="PCA", ncomponents=4, colour_by="total_features") + fontsize + ggtitle("Reduce Dimensions with PCA")

out1 <- plotTSNE(sce, use_dimred="PCA", perplexity=4, colour_by="celltypes", 
                 rand_seed=100) + fontsize + ggtitle("Perplexity = 4")
out2 <- plotTSNE(sce, use_dimred="PCA", perplexity=5, colour_by="celltypes",
                  rand_seed=100) + fontsize + ggtitle("Perplexity = 5")
out3 <- plotTSNE(sce, use_dimred="PCA", perplexity=6, colour_by="celltypes",
                  rand_seed=100) + fontsize + ggtitle("Perplexity = 6")
out4 <- plotTSNE(sce, use_dimred="PCA", perplexity=7, colour_by="celltypes",
                 rand_seed=100) + fontsize + ggtitle("Perplexity = 7")

multiplot(out1, out2, out3, out4, cols=2)

Brennecke_HVG <- BrenneckeGetVariableGenes(
  expr_mat = normcounts(sce),
  spikes = NA,
  fdr = 0.1,
  minBiolDisp = 0.25
)

## another method to identify by Kolodziejczyk AA, Kim JK, Tsang JCH et al. (2015)
means <- rowMeans(normcounts(sce))
cv2 <- apply(normcounts(sce), 1, var)/means^2
dm.stat <- DM(means, cv2)
#head(dm.stat)


DM_HVG = names(dm.stat)[which(dm.stat>0.3)]
#DM_HVG = DM_HVG[which()]
dev.off()

sce.HVG.Brenneck = sce[rownames(sce)%in%Brennecke_HVG, ]
sce.HVG.DM = sce[rownames(sce)%in%DM_HVG, ]

save(sce, sce.HVG.Brenneck, sce.HVG.DM, file=paste0('../Rdata/', version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
#HVG_genes <- Brennecke_HVG$Gene

##################################################
## clustering parts (test different methods)
##################################################
library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
set.seed(1234567)
library(dynamicTreeCut)

load(file=paste0('../Rdata/', version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))

### test workflow by Aaron T. L. Lun 
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))

#sce.sels = sce.HVG.Brenneck
sce.sels = sce.HVG.DM
sce.sels = sce.sels[, which(sce.sels$celltypes != "ASEL.bulk" & sce$celltypes != "ASER.bulk")]

plotTSNE(sce.sels, colour_by="celltypes", size_by = "total_features",
         perplexity=5, rand_seed=100)

 
TEST.Aaron.workflow = TRUE
  if(TEST.Aaron.workflow)
  {  
  #pcs <- reducedDim(sce, "PCA")
  pcs = t(logcounts(sce.sels))
  my.dist <- dist(pcs)
  my.tree <- hclust(my.dist, method="ward.D2")
  
  my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0, minClusterSize = 2))
  sce.sels$cluster <- factor(my.clusters)
 
  plotTSNE(sce.sels, colour_by="cluster", shape_by = "celltypes", size_by = "total_features",
           perplexity=6, rand_seed=100) + fontsize
  
}

### test SC3
TEST.SC3 = TRUE
if(TEST.SC3)
{
  sce = sce.sels
  # create a SingleCellExperiment object
  sce <- SingleCellExperiment(
    assays = list(
      counts = counts(sce),
      logcounts = log2(normcounts(sce) + 1)
    ), 
    colData = colData(sce),
    rowData = rowData(sce)
  )
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  # remove features with duplicated names
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  
  # define spike-ins
  #isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)
  
  
  #var.fit <- trendVar(sce, parametric=TRUE, span=0.2)
  plotPCA(sce.sels,
          exprs_values = "logcounts",
          colour_by = "celltypes",
          size_by = "total_features"
  )
  
  ### test SC3 clustering method
  sce <- sc3(sce.sels, gene_filter = FALSE, ks = 2:6, biology = TRUE)
  
  #sce = sc3_estimate_k(sce)
  #metadata(deng)$sc3$k_estimation
  #rowData(sce)$feature_symbol
  #sce <- sc3(sce, ks = 2, gene_filter = TRUE, biology = TRUE)
  #deng <- plotTSNE(deng, rand_seed = 1, return_SCE = TRUE)
  
  col_data <- colData(sce)
  head(col_data[ , grep("sc3_", colnames(col_data))])
  plotPCA(
    sce, 
    colour_by = "sc3_3_clusters", 
    size_by = "sc3_3_log2_outlier_score"
  )
  
  plotPCA(
    sce, 
    colour_by = "sc3_3_clusters",
    shape_by = "celltypes",
    size_by = "total_features"
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
  
  plotTSNE(sce, colour_by="sc3_3_clusters", shape_by = "celltypes", size_by = "total_features",
           perplexity=6, rand_seed=100) + fontsize
  
  #sc3_plot_consensus(sce, k = 3)
  sc3_plot_consensus(
    sce, k = 3, 
    show_pdata = c(
      "celltypes", 
      "sc3_3_clusters", 
      "log10_total_features",
      "log10_total_counts",
      
      "sc3_3_log2_outlier_score"
    )
  )
  
  sc3_plot_cluster_stability(sce, k = 6)
  
  sc3_plot_de_genes(sce, k = 3, p.val = 0.2)
  
  sc3_plot_markers(sce, k = 3)
  
}


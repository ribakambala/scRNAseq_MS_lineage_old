##################################################
##################################################
## Project: transcriptional priming mechanism by tbx in C.elegans
## Script purpose: analysis the single-cell RNA-seq data
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Feb 19 14:43:38 2018
##################################################
##################################################
version.DATA = 'R6875_R7116_R7130_scRNA_v1'
version.analysis =  paste0(version.DATA, '_20190116')

design.file = "../exp_design/R6875_sample_infos.xlsx"
dataDir = paste0("../data/")

resDir = paste0("../results/", version.analysis)
tabDir = paste0("../results/", version.analysis, "/tables/")
RdataDir = paste0("../results/", version.analysis, "/Rdata/")
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

Manually.make.design.matrix = TRUE
Use.sampleID.mapSamples = FALSE
##################################################
##################################################
## Section: Import data first 
## and then add metadata from the design matrix 
## and then the quality controls from nf-RNAseq
##################################################
##################################################
library(data.table)
# first concatenate the data
xlist <-list.files(path=dataDir, pattern = "*merged_gene_counts.txt", full.names = TRUE)
aa = NULL

for(n in 1:length(xlist))
{
  # n = 1
  cat(n, '\t')
  cat(xlist[n], '\n')
  
  if(n==1){
    #aa = read.delim(xlist[n], header = TRUE);
    #colnames(aa) = c('RefseqID', basename(xlist[n]))
    #aa = data.frame(aa, stringsAsFactors = FALSE)
    aa = data.frame(fread(xlist[n], header=TRUE, sep="\t", stringsAsFactors=FALSE), stringsAsFactors = FALSE)
  }else{
    test = data.frame(fread(xlist[n], header=TRUE, sep="\t", stringsAsFactors=FALSE), stringsAsFactors = FALSE)
    mm = match(aa$ENSEMBL_ID, test$ENSEMBL_ID)
    aa = data.frame(aa, test[mm, -grep('ENSEMBL_ID', colnames(test))])
  }
}

colnames(aa)[1] = 'gene';
jj = grep("CCVTBANXX_8.76090_", colnames(aa))
aa = aa[, c(1, jj)]
colnames(aa)[-1] = sapply(colnames(aa)[-1], function(x) gsub("CCVTBANXX_8.76090_", "", x))

if(Use.sampleID.mapSamples){
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
}else{
  
  xx = colnames(aa)[-1]
  xx = data.frame(xx, rep("others", length(xx)), stringsAsFactors = FALSE)
  colnames(xx) = c("barcodes", "sampleInfo")  
  mm = match(xx$barcodes, design$barcodes)
  xx$sampleInfo[which(!is.na(mm))] = design$sampleInfo[mm[which(!is.na(mm))]]
  design = xx;
  
}

library("openxlsx")
design = read.xlsx(design.file, sheet = 1, colNames = TRUE)

if(Manually.make.design.matrix){
  library(stringi)
  design = data.frame(design$`Multiplex/Barcode`, design$Brief.Sample.Description, stringsAsFactors = FALSE)
  colnames(design) = c("barcodes", "sampleInfo")
  xx = design
  design = xx
  design$barcodes = gsub('[[:digit:]]+', '', design$barcodes)
  design$barcodes = gsub(':', '', design$barcodes, fixed = TRUE)
  design$barcodes = gsub("[ +]", "", design$barcodes)
  design$barcodes = gsub("\\s", "", design$barcodes)    
  design = design[-c(1:2), ]
  design$sampleInfo[grep("bulk", design$sampleInfo)] = "bulk.control"
  design$sampleInfo[grep("early", design$sampleInfo)] = "early"
  
  #design$barcodes = stri_replace_all(string=design$barcodes, pattern=" ", repl="", fixed = TRUE)
   
    #design$barcodes = gsub("[^*:]", "", design$barcodes)
  #design$barcodes = gsub(" ", "", design$barcodes)
  #design = design[which(!is.na(design$sampleID)), ]
  #design$celltypes = sapply(design$celltypes, function(x) gsub("\\(", "", as.character(x)))
  #design$celltypes = sapply(design$celltypes, function(x) gsub("\\)", "", as.character(x)))
  #design$celltypes = sapply(design$celltypes, function(x) gsub(" ", ".", as.character(x)))
  
  #conditions = sapply(design$sampleInfo, function(x) paste0(rev(rev(unlist(strsplit(as.character(x), "_")))[-c(1:3)]), collapse="_"))
  #times =  sapply(design$sampleInfo, function(x) rev(unlist(strsplit(as.character(x), "_")))[3])
  
  #design = data.frame(design[, c(1, 3, 2)],stringsAsFactors = FALSE)
  #design$conditions = sapply(design$conditions, function(x) gsub("_", "-", x))
}



save(design, aa, file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_RNA_seq.Rdata'))

##################################################
##################################################
# Section: Quality control, process and clean the count table for scRNA-seq
# 1) remove genes with zero read across all cells
# 2) convert emsemble gene annotation to gene symbols
# 3) make singleCellExperiment object
##################################################
##################################################
load( file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_RNA_seq.Rdata'))
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
#library(limma)
options(stringsAsFactors = FALSE)

## make SCE object and remove genes with zero reads detected
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = as.data.frame(design), 
                            rowData = data.frame(gene_names = rownames(counts), feature_symbol = rownames(counts)))

write.csv(counts(sce), file=paste0(tabDir, "scRNAseq_raw_readCounts", version.analysis, ".csv"), row.names=TRUE)

keep_feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_feature, ]

#is.spike <- grepl("^ERCC", rownames(sce))
is.mito <- rownames(sce) %in% gg.Mt;
is.ribo <- rownames(sce) %in% gg.ribo;
summary(is.mito)

sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito, Ribo=is.ribo))
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

hist(sce$log10_total_counts, xlab="Library sizes (in log10)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
abline(v=c(6, 5, 4), col="red", lwd=2.0)

hist(sce$total_features_by_counts_endogenous, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
#hist(sce$pct_counts_ERCC, xlab="ERCC proportion (%)", 
#     ylab="Number of cells", breaks=20, main="", col="grey80")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

hist(sce$pct_counts_Ribo, xlab="ribo proportion (%)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

par(mfrow=c(1, 1))
plot(sce$total_counts, sce$total_features_by_counts, xlab="Library size", ylab='Nb of expressed genes', log='xy')
abline(v=c(10^4, 10^5, 10^6), col='red', lwd=2.0)
abline(h=c(100, 200, 1000, 2000), col='red', lwd=2.0)

plotColData(sce, 
            y = "log10_total_features_by_counts",
            x = "log10_total_counts",
            colour_by = "pct_counts_Ribo",
            size_by = "pct_counts_Mt"
)

dev.off()
## here we are using the 50,000 for library size and 100 expressed genes as thresholds
threshod.total.counts.per.cell = 10^5
threshod.nb.detected.genes.per.cell = 1000;
#libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
#feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
#libsize.drop = sce$total_counts<5*10^4
#feature.drop = sce$total_features < 200
filter_by_total_counts <- (sce$total_counts > threshod.total.counts.per.cell)
table(filter_by_total_counts)
filter_by_expr_features <- (sce$total_features_by_counts > threshod.nb.detected.genes.per.cell)
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
#sce = normalize(sce)

#fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
#plotPCA(sce, pca_data_input="pdata") + fontsize

## comparison between manual filtering and automatic ouliter filtering
Manual.vs.outlier.filtering = FALSE
if(Manual.vs.outlier.filtering)
{
  sce <- plotPCA(
    sce,
    size_by = "total_features_by_counts", 
    #shape_by = "",
    by_exprs_values = "logcounts",
    shape_by = "use",
    return_SCE = FALSE
  )
  
  table(sce$outlier)
  
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

sce = sce[, sce$use]
save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_filtered_SCE.Rdata'))

####################
## filter lowly expressed (and probably too highly expressed genes)
####################
load(file=paste0(RdataDir, version.DATA, '_QCed_cells_filtered_SCE.Rdata'))

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

genes.to.keep <- num.cells > 5 & ave.counts >= 1  & ave.counts <10^6  # detected in >= 2 cells, ave.counts >=5 but not too high
summary(genes.to.keep)

# remove mt and ribo genes
genes.to.keep = genes.to.keep & ! rownames(sce) %in% gg.Mt & ! rownames(sce) %in% gg.ribo
summary(genes.to.keep)

sce <- sce[genes.to.keep, ]

save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_SCE.Rdata'))

########################################################
########################################################
# Section : chekc confonding factors, normalization and 
# batch correction
# codes original from Hemberg's course
# https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#data-visualization
##########
# Among all normalization methods, DESeq2 and scran normalization will be used
# consider seurat for UMI counts, but here we are just working on the read counts
########################################################
########################################################
Data.visulization.exploration = FALSE
if(Data.visulization.exploration)
{
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_SCE.Rdata'))
  sce.qc = sce
  
  assay(sce, "logcounts_raw") <- log2(counts(sce) + 1)
  reducedDim(sce) <- NULL
  #reducedDim(sce.qc) <- NULL
  endog_genes <- !rowData(sce)$is_feature_control
  
  pdfname = paste0(resDir, "/scRNAseq_filtered_identifying_confounding_factors_", version.analysis, ".pdf")
  pdf(pdfname, width=10, height = 6)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  library(scater, quietly = TRUE)
  options(stringsAsFactors = FALSE)
  
  plotPCA(
    sce[endog_genes, ],
    run_args = list( exprs_values = "logcounts_raw"), 
    colour_by = "total_counts",
    size_by = "total_features_by_counts"
    #shape_by = "individual"
  )
  
  plotDiffusionMap(
    sce[endog_genes, ],
    run_args = list( exprs_values = "logcounts_raw"), 
    colour_by = "total_counts",
    size_by = "total_features_by_counts"
    #shape_by = "individual"
  )
  
  scater::plotMDS(
    sce[endog_genes, ],
    run_args = list( exprs_values = "logcounts_raw"), 
    colour_by = "total_counts",
    size_by = "total_features_by_counts"
    #shape_by = "individual"
  )
  
  plotUMAP(
    sce[endog_genes, ],
    run_args = list( exprs_values = "logcounts_raw"), 
    colour_by = "total_counts",
    size_by = "total_features_by_counts"
    #shape_by = "individual"
  )
 
  param.perplexity = 6;
  set.seed(1)
  plotTSNE(
    sce, 
    run_args=list(exprs_values="logcounts_raw", 
                  perplexity=param.perplexity),
    colour_by = "total_counts",
    size_by = "total_features_by_counts"               
  )
  
  
  ##### identify confounding factors
  plotQC(
    sce[endog_genes, ],
    type = "find-pcs",
    exprs_values = "logcounts_raw",
    variable = "total_features_by_counts"
  )
  
  plotQC(
    sce[endog_genes, ],
    type = "expl",
    exprs_values = "logcounts_raw",
    variables = c(
      "total_features_by_counts",
      "total_counts",
      "pct_counts_Mt"
    )
  )
  
  dev.off()

}

##################################
# scRNA-seq data normalization
# here many normalization have been proposed
# but here we are still using TMM from edgeR or from scran
##################################
load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_SCE.Rdata'))
assay(sce, "logcounts_raw") <- log2(counts(sce) + 1)
reducedDim(sce) <- NULL
sce.qc = sce
endog_genes <- !rowData(sce)$is_feature_control

# specify the shapes for bulk control and early cells (supposed to be)
shapes = rep(1, ncol(sce))
shapes[which(sce$sampleInfo=="bulk.control")] = 3
shapes[which(sce$sampleInfo=="early")] = 2
shapes = data.frame(shapes)

library(scran)
options(stringsAsFactors = FALSE)
set.seed(1234567)

pdfname = paste0(resDir, "/scRNAseq_filtered_normalization_", version.analysis, ".pdf")
pdf(pdfname, width=10, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

### raw counts in log scale
plotPCA(
  sce.qc[endog_genes, ],
  run_args = list(exprs_values = "logcounts_raw"), 
  colour_by = "sampleInfo",
  size_by = "total_features_by_counts"
  #shape_by = "sampleInfo"
) + ggtitle("PCA using logcounts_raw")

### cpm
logcounts(sce.qc) <- log2(calculateCPM(sce.qc, use_size_factors = FALSE) + 1)
scater::plotPCA(
  sce.qc[endog_genes, ],
  run_args = list(exprs_values = "logcounts"), 
  colour_by = "total_counts",
  size_by = "total_features_by_counts"
  #shape_by = "individual"
) + ggtitle("PCA using cpm")

try.TMM = FALSE
if(try.TMM){
  ### TMM
  sce.qc <- normaliseExprs(
    sce.qc,
    method = "TMM",
    feature_set = endog_genes,
    return_log = TRUE,
    return_norm_as_exprs = TRUE
  )
  plotPCA(
    sce.qc[endog_genes, ],
    colour_by = "total_counts",
    size_by = "total_features_by_counts"
  )
}

## scran normalization (not working here, because negative scaling factor found)
qclust <- quickCluster(sce.qc, min.size = 5)
sce.qc <- computeSumFactors(sce.qc, sizes = 5, clusters = qclust)
sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)

plotPCA(
  sce.qc[endog_genes, ],
  colour_by = "total_counts",
  size_by = "total_features_by_counts",
  shape_by = shapes
) + ggtitle("scran normalization")

param.perplexity = 10;
plotTSNE(
  sce.qc[endog_genes, ],
  run_args = list(exprs_values = "logcounts", perplexity = param.perplexity), 
  colour_by = "total_counts",
  size_by = "total_features_by_counts",
  shape_by = shapes
) + ggtitle(paste0("tSNE - perplexity = ", param.perplexity))

if(any(sizeFactors(sce.qc)<0)){
  cat("ERROR........")
  cat("NEGATIVE size factors FOUND !!! \n" )
  cat(sizeFactors(sce), "\n")
}else{
  plot(sizeFactors(sce.qc), sce.qc$total_counts/1e6, log="xy",
       ylab="Library size (millions)", xlab="Size factor")
  
}

dev.off()

## save the normalized sce object
sce = sce.qc;
save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))

write.csv(logcounts(sce), file=paste0(tabDir, "scRNAseq_QC_filtered_normalized_readCounts.csv"), row.names=TRUE)

########################################################
########################################################
# Section : 
# feature selection and dimension reduction and clustering
# there are many options for HGV, dimension reduction and clustering
# there are also many choices to combine them
# Here, to start with something simple, I test three popular workflow from Aaron, hemberg and seurat
########################################################
########################################################
library(scater)
library(scran)
library(scRNA.seq.funcs)
library(matrixStats)
library(M3Drop)
library(RColorBrewer)
library(SingleCellExperiment)
set.seed(1)

load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))

##########################################
# test Aaron's clustering workflow
# original code found in 
# https://master.bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/work-1-reads.html/
# #7_modelling_the_technical_noise_in_gene_expression
##########################################
TEST.Aaron.scRNAseq.workflow.clustering.part = FALSE
if(TEST.Aaron.workflow)
{ 
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
  
  var.fit <- trendVar(sce,method="loess", parametric=TRUE, use.spikes = FALSE, assay.type="logcounts",
                      loess.args=list(span=0.3))
  var.out <- decomposeVar(sce, var.fit)
  head(var.out)
  
  #plot(var.fit.nospike$mean, var.fit.nospike$var)
  
  pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_HVG_featureSelect_clustering_", version.analysis, ".pdf")
  pdf(pdfname, width=10, height = 6)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  par(mfrow=c(1, 1))
  plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
       ylab="Variance of log-expression")
  curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
  #cur.spike <- isSpike(sce)
  #points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
  
  chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:50]
  plotExpression(sce, features=rownames(var.out)[chosen.genes]) + fontsize
  #chosen.genes <- order(var.out.nospike$bio, decreasing=TRUE)[1:50]
  #plotExpression(sce, features=rownames(var.out.nospike)[chosen.genes]) + fontsize
  
  sce <- denoisePCA(sce, technical=var.out)
  dim(reducedDim(sce, "PCA")) 
  
  #sce <- denoisePCA(sce, technical=var.fit.nospike$trend) 
  #dim(reducedDim(sce, "PCA")) 
  
  # low-dimension visulization
  plotReducedDim(sce, use_dimred="PCA", ncomponents=4, 
                 colour_by="total_features_by_counts") + fontsize
  
  #plotReducedDim(sce, use_dimred="PCA", ncomponents=4, colour_by="total_features_by_counts") + fontsize + 
  
  set.seed(100)
  out5 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=5),
                   colour_by="total_features_by_counts") + fontsize + ggtitle("Perplexity = 5")
  
  set.seed(100)
  out10 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=10),
                    colour_by="total_features_by_counts") + fontsize + ggtitle("Perplexity = 10")
  
  set.seed(100)
  out20 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=20),
                    colour_by="total_features_by_counts") + fontsize + ggtitle("Perplexity = 20")
  
  multiplot(out5, out10, out20, cols=3)
  
  set.seed(100)
  sce <- runTSNE(sce, use_dimred="PCA", perplexity=10)
  reducedDimNames(sce)
  
  ##########################################
  # clustering part 
  ##########################################
  pcs <- reducedDim(sce, "PCA")
  my.dist <- dist(pcs)
  my.tree <- hclust(my.dist, method="ward.D2")
  
  hist(my.tree)
  
  library(dynamicTreeCut)
  my.clusters <- unname(cutreeDynamic(my.tree, method = "hybrid", 
                                      distM=as.matrix(my.dist), 
                                      minClusterSize=2, verbose=1))
  
  sce$cluster <- factor(my.clusters)
  set.seed(100)
  plotTSNE(sce,
           run_args = list(perplexity = 20),
           colour_by="cluster", size_by = "total_features_by_counts", shape_by = shapes) + fontsize
  
  library(cluster)
  clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
  sil <- silhouette(my.clusters, dist = my.dist)
  sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
  sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
  plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
       border=sil.cols, col=sil.cols, do.col.sort=FALSE) 
  

  markers <- findMarkers(sce, my.clusters)
  
  marker.set <- markers[["1"]]
  head(marker.set, 10)
  
  top.markers <- rownames(marker.set)[marker.set$Top <= 5]
  plotHeatmap(sce, features=top.markers, 
              columns=order(sce$cluster), 
              colour_columns_by=c("cluster"),
              cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5)) 
  
  dev.off()
  
}

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


########################################################
########################################################
# Section : pseudo-time analysis
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







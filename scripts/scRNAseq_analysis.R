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


dataDir = paste0("../data/")

resDir = paste0("../results/", version.analysis)
tabDir = paste0("../results/", version.analysis, "/tables/")
RdataDir = paste0("../results/", version.analysis, "/Rdata/")
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

Manually.Specify.sampleInfos.4scRNAseq = TRUE
Aggregate.nf.QCs.plots.in.designMatrix = TRUE
#design.file = "../exp_design/R6875_sample_infos.xlsx"
#Use.sampleID.mapSamples = FALSE
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
    aa = data.frame(fread(xlist[n], header=TRUE, sep="\t", stringsAsFactors=FALSE), stringsAsFactors = FALSE)
  }else{
    test = data.frame(fread(xlist[n], header=TRUE, sep="\t", stringsAsFactors=FALSE), stringsAsFactors = FALSE)
    mm = match(aa$ENSEMBL_ID, test$ENSEMBL_ID)
    aa = data.frame(aa, test[mm, -grep('ENSEMBL_ID', colnames(test))])
  }
}

colnames(aa)[1] = 'gene';

##########################################
# manually make design matrix 
##########################################
if(Manually.Specify.sampleInfos.filtering.4scRNAseq){
  library(stringi)
  library("openxlsx")
  
  design = data.frame(samples = colnames(aa)[-1], stringsAsFactors = FALSE)
  design$flowcell.lane = gsub(".\\w+$", "",  design$samples)
  design$sampleIDs = gsub('^\\w+_\\d.', "", design$samples)
  design$sampleIDs = gsub('_\\w+$', '', design$sampleIDs)
  design$barcodes = gsub('^\\w+_\\d.\\d+_', '', design$samples)
  
  ##########################################
  # here manually 
  ##########################################
  design$request = NA;
  design$request[which(design$flowcell.lane == "CCVTBANXX_8")] = "R6875"
  design$request[which(design$flowcell.lane == "CCVBPANXX_1")] = "R7116"
  design$request[which(design$flowcell.lane == "HHG5KBGX9_1")] = "R7130"
  design$request[which(design$flowcell.lane == "HHGHNBGX9_1")] = "R7130"
  design$seqInfos = paste0(design$request, "_", design$flowcell.lane)
  #jj = grep("CCVTBANXX_8.76090_", colnames(aa))
  #aa = aa[, c(1, jj)]
  #colnames(aa)[-1] = sapply(colnames(aa)[-1], function(x) gsub("CCVTBANXX_8.76090_", "", x))
  kk = which(!is.na(design$request))
  design = design[kk, ]
  aa = aa[, c(1, kk+1)]
  
}else{
  #design = read.xlsx(design.file, sheet = 1, colNames = TRUE)
  #design = data.frame(design$`Multiplex/Barcode`, design$Brief.Sample.Description, stringsAsFactors = FALSE)
  #colnames(design) = c("barcodes", "sampleInfo")
  #xx = design
  # design = xx
  # design$barcodes = gsub('[[:digit:]]+', '', design$barcodes)
  # design$barcodes = gsub(':', '', design$barcodes, fixed = TRUE)
  # design$barcodes = gsub("[ +]", "", design$barcodes)
  # design$barcodes = gsub("\\s", "", design$barcodes)    
  # design = design[-c(1:2), ]
  # design$sampleInfo[grep("bulk", design$sampleInfo)] = "bulk.control"
  # design$sampleInfo[grep("early", design$sampleInfo)] = "early"
  # 
  # if(Use.sampleID.mapSamples){
  #   index = c()
  #   for(n in 1:nrow(design))
  #   {
  #     jj = grep(design$sampleID[n], colnames(aa))
  #     if(length(jj)==1){index = c(index, jj)
  #     }else{cat("NOT FOUND sample for ", design$sampleID[n], "\n")}
  #   }
  #   
  #   aa = data.frame(aa$gene, aa[, index], stringsAsFactors = FALSE)
  #   colnames(aa)[-1] = apply(design[, c(2,1)], 1, paste0, collapse="_")
  #   colnames(aa)[1] = 'gene';
  # }else{
  #   
  #   xx = colnames(aa)[-1]
  #   xx = data.frame(xx, rep("others", length(xx)), stringsAsFactors = FALSE)
  #   colnames(xx) = c("barcodes", "sampleInfo")  
  #   mm = match(xx$barcodes, design$barcodes)
  #   xx$sampleInfo[which(!is.na(mm))] = design$sampleInfo[mm[which(!is.na(mm))]]
  #   design = xx;
  #   
  # }
}

##########################################
# aggregated quality controls from nf-RNAseq 
##########################################
if(Aggregate.nf.QCs.plots.in.designMatrix){
  #load(file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_RNA_seq.Rdata'))
  
  source('functions_aggregate_nf_qc.R')
  dirs.all = c('../../../Ariane/R7116_R7130_scrnaseq/results_all/multiqc_data_1', 
               "../../../Ariane/R6875_scRNAseq/results_all_3rd/MultiQC/multiqc_data")
  QCs.nf = aggrate.nf.QCs(dirs.all)
  QCs.nf$Sample = gsub("#", ".", QCs.nf$Sample)
  
  mm = match(design$samples, QCs.nf$Sample)
  xx = data.frame(design, QCs.nf[mm, ], stringsAsFactors = FALSE)
  
  design = xx;
   
}

save(aa, design, file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_design_sampleInfos_QCs_nf_RNA_seq.Rdata'))

##################################################
##################################################
# Section: Quality control, process and clean the count table for scRNA-seq
##################################################
##################################################
load(file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_design_sampleInfos_QCs_nf_RNA_seq.Rdata'))

##########################################
# first round cleaning the count table:
# a) keep only rows for mapped transripts
# b) map the gene names to the gene symbols
##########################################
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

##########################################
# Import SingleCellExperiment and scater packages for the QC and table cleaning
# several steps will be proceded:
# 1) general overview of data quality: sequencing depth, mapping rate, assignment rate, rRAN codamination for each sequencing lane
# 2) clean the cells 
# 3) clean genes 
##########################################
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)

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

####################
## filter cells with low quality 
####################
pdfname = paste0(resDir, "/scRNAseq_QCs_cells_filterting.pdf")
pdf(pdfname, width=10, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

# some general statistics for each request and lane
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

dev.off()

##########################################
## filter cells with low quality 
# here we are using the 50,000 for library size and 100 expressed genes as thresholds
##########################################
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
filter_by_MT = sce$pct_counts_Mt < 7.5
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
if(Manual.vs.outlier.filtering){
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

##########################################
##########################################
# remove the cell cycle confounder 
# we need to train the cells to identify the cell cycle phase
# this could be more complicated than expected
##########################################
##########################################
set.seed(100)
library(scran)

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                package="scran"))

assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
plot(assignments$score$G1, assignments$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16)

sce$phases <- assignments$phases
table(sce$phases)

####################
## filter lowly expressed (and probably too highly expressed genes)
####################
load(file=paste0(RdataDir, version.DATA, '_QCed_cells_filtered_SCE.Rdata'))

pdfname = paste0(resDir, "/scRNAseq_QCs_genes_filterting_", version.analysis, ".pdf")
pdf(pdfname, width=10, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

#plotQC(reads, type = "highest-expression", n=20)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotHighestExprs(sce, n=50) + fontsize

#fontsize <- theme(axis.text=element_text(size=16), axis.title=element_text(size=16))
#plotHighestExprs(sce, n=30) + fontsize

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
# Section : normalization and confonding factors,
# batch correction
# codes original from Hemberg's course
# https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#data-visualization
##########
# Among all normalization methods, DESeq2 and scran normalization will be used
# consider seurat for UMI counts, but here we are just working on the read counts
########################################################
########################################################
##################################
# scRNA-seq data normalization 
# Many normalization have been proposed
# Here we test main two methods: TMM from edgeR or from DESeq2
# and the one from scran
# the PCA and some other plots were used to assess the normalization
##################################
load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_SCE.Rdata'))
library(scRNA.seq.funcs)
library(scater)
library(scran)
options(stringsAsFactors = FALSE)

sce.qc = sce;
reducedDim(sce) <- NULL
endog_genes <- !rowData(sce)$is_feature_control

Methods.Normalization = c("cpm", "DESeq2", "UQ", "scran",
                          "downsample")

#Methods.Normalization = "DESeq2" 

for(method in Methods.Normalization)
{
  
  # method = Methods.Normalization[4]
  set.seed(1234567)
  
  cat('normalization method -- ', method, "\n")
  
  pdfname = paste0(resDir, "/scRNAseq_filtered_normalization_testing_", method, ".pdf")
  pdf(pdfname, width=14, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  main = paste0(method, " normalization");
  if(method == "raw") { # raw log counts
    assay(sce.qc, "logcounts") <- log2(counts(sce.qc) + 1)
  }
  
  cat("start to normalize the data ---\n")
  
  if(method == "cpm") { ### cpm
    assay(sce.qc, "logcounts") <- log2(calculateCPM(sce.qc, use_size_factors = FALSE) + 1)
  }
  if(method == "UQ"){
    source("scRNAseq_functions.R")
    logcounts(sce.qc) <- log2(cal_uq_Hemberg(counts(sce.qc)) + 1)
  }
  
  if(method == "DESeq2"){
    source("scRNAseq_functions.R")
    sizeFactors(sce.qc) = calculate.sizeFactors.DESeq2(counts(sce.qc))
    sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
  }
  if(method == "downsample") {
    assay(sce.qc, "logcounts") <- log2(Down_Sample_Matrix(counts(sce.qc)) + 1)
  }
  
  if(method == "scran"){
    ## scran normalization (not working here, because negative scaling factor found)
    qclust <- quickCluster(sce.qc, min.size = 30)
    sce.qc <- computeSumFactors(sce.qc, sizes = 15, clusters = qclust)
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
  
  scater::plotPCA(
    sce.qc[endog_genes, ],
    run_args = list(exprs_values = "logcounts"), 
    size_by = "total_counts",
    #size_by = "total_features_by_counts",
    colour_by = "seqInfos"
  ) + ggtitle(paste0("PCA -- ", main))
  
  plotRLE(
    sce.qc[endog_genes, ],
    exprs_values = "logcounts", 
    exprs_logged = TRUE,
    colour_by = "seqInfos"
  ) + ggtitle(paste0("RLE -- ", main))
  
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
      
  #plotDiffusionMap(
  #  sce.qc[endog_genes, ],
  #  run_args = list(exprs_values = "logcounts"), 
  #  size_by = "total_counts",
    #size_by = "total_features_by_counts",
  #  colour_by = "seqInfos"
  #) + ggtitle(paste0("DM -- ", main))
  
  dev.off()
  
  
}


## select normalization method and save the normalized sce object
set.seed(1000)    
#qclust <- quickCluster(sce.qc, min.size = 30)
clusters <- quickCluster(sce, method="igraph", min.mean=0.1)
table(clusters)

sce <- computeSumFactors(sce, min.mean=0.1, clusters = clusters)
sce <- normalize(sce, exprs_values = "counts", return_log = TRUE)

#sce = sce.qc;
save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))


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
load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))

library(scRNA.seq.funcs)
#library(RUVSeq)
library(scater)
library(SingleCellExperiment)
library(scran)
#library(kBET)
#library(sva) # Combat
#library(edgeR)
set.seed(1234567)
options(stringsAsFactors = FALSE)

#sce.qc = sce

##########################################
# here to prepare the inputs for MNN in scran 
# Feature selection (HVGs) is performed
# here we just use the scran trandVar()
# further test to follow for other methods
# (https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby011/4898116)
# But many of them, e.g. BASics, Brennecke works better with spike-in 
##########################################
#block <- paste0(sce.gse81076$Plate, "_", sce.gse81076$Donor)
#fit <- trendVar(sce.gse81076, block=block, parametric=TRUE)
#dec <- decomposeVar(sce.gse81076, fit)
block <- sce$request
fit <- trendVar(sce, block=block, parametric=TRUE, assay.type="logcounts", use.spikes=FALSE)
dec <- decomposeVar(sce, fit)

plot(dec$mean, dec$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)

dec.sorted <- dec[order(dec$bio, decreasing=TRUE),]
head(dec.sorted)
#length(which(dec.sorted$bio>0))

# here HVGs selected with FDR<0.01
gene.chosen <- rownames(dec.sorted)[which(dec.sorted$FDR <0.01)]
length(gene.chosen)

##########################################
# Here we are going to run both mnnCorrect() and fastMNN() 
##########################################
rescaling.for.multBatch = FALSE
if(rescaling.for.multBatch){
  rescaled <- multiBatchNorm(sce[ , which(sce.qc$request == "R6875")], 
                             sce[ , which(sce.qc$request == "R7116")],
                             sce[ , which(sce.qc$request == "R7130")])
  rescaled.R6875 <- rescaled[[1]]
  rescaled.R7116 <- rescaled[[2]]
  rescaled.R7130 <- rescaled[[3]] 
  
  R6875=logcounts(rescaled.R6875)
  R7116=logcounts(rescaled.R7116)
  R7130=logcounts(rescaled.R7130)
}else{
  R6879 = logcounts(sce[ , which(sce$request == "R6875")])
  R7116 = logcounts(sce[ , which(sce$request == "R7116")])
  R7130 = logcounts(sce[ , which(sce$request == "R7130")])
}

Use.fastMNN = TRUE
if(Use.fastMNN){
  set.seed(100) 
  original <- list(
    R6879 = R6879,
    R7116 = R7116,
    R7130 = R7130
  )
  
  # Slightly convoluted call to avoid re-writing code later.
  # Equivalent to fastMNN(GSE81076, GSE85241, k=20, d=50, approximate=TRUE)
  mnn.out <- do.call(fastMNN, c(original, list(k=20, cos.norm = TRUE, d=50, subset.row = gene.chosen, auto.order=c(3,2,1),
                                               approximate=TRUE)))
  dim(mnn.out$corrected)
  mnn.out$batch
  mnn.out$pairs
  
  #omat <- do.call(cbind, original)
  #sce.qc <- SingleCellExperiment(list(logcounts=omat))
  reducedDim(sce, "MNN") <- mnn.out$corrected
  sce$Batch <- as.character(mnn.out$batch)
  sce
  
  set.seed(1000)
  with.var <- do.call(fastMNN, c(original,
                                 list(k=20, cos.norm = TRUE, d=50, subset.row = gene.chosen, auto.order=c(3,2,1), approximate=TRUE,
                                      compute.variances=TRUE)))
  with.var$lost.var
  
  #save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_batchCorrect_fastMNN.Rdata'))
  
  pdfname = paste0(resDir, "/scRNAseq_filtered_test_batchCorrection_inLowDimensions.pdf")
  pdf(pdfname, width=18, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  # check the effect of correction by fastMNN
  set.seed(100)
  # Using irlba to set up the t-SNE, for speed.
  osce <- runPCA(sce, ncomponents = 50, ntop=Inf, method="irlba", exprs_values = "logcounts", feature_set = gene.chosen)
  #plotPCA(osce, colour_by="Batch") + ggtitle("Original")
  osce <- runTSNE(osce, use_dimred="PCA", perplexity = 20)
  ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original")
  set.seed(100)
  csce <- runTSNE(sce, use_dimred="MNN", perplexity = 20)
  ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")
  multiplot(ot, ct, cols=2)
  
  osce <- runUMAP(osce, use_dimred="PCA", perplexity = 20)
  ot <- plotUMAP(osce, colour_by="Batch") + ggtitle("Original")
  set.seed(100)
  csce <- runUMAP(sce, use_dimred="MNN", perplexity = 20)
  ct <- plotUMAP(csce, colour_by="Batch") + ggtitle("Corrected")
  multiplot(ot, ct, cols=2)
  
  dev.off()
  
}

Use.mnnCorrect = TRUE
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
  
  pbmc = FindClusters(pbmc, resolution = 4, algorithm = 3)
  sce$cluster <- factor(pbmc@active.ident)
  
  plotTSNE(sce, colour_by="cluster", size_by = "total_features_by_counts") + fontsize + ggtitle("seurat - graph base clustering")
  plotUMAP(sce, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
    fontsize + ggtitle("seurat -- graph based clustering")
  
  Select.Early.timePoints = FALSE
  if(Select.Early.timePoints){
     
    xx = table(sce$cluster,sce$Batch)
    #colnames(xx) = sce$Batch
    cluster4early = rownames(xx)[which(xx[, 1]>=5|xx[,2]>=5)]
    
    mm = match(sce$cluster, factor(cluster4early))
    
    sels = which(!is.na(mm))
    
    sce.sel = sce[, sels ]
    set.seed(100)
    sce.sel <- runTSNE(sce.sel,  use_dimred = "MNN", perplexity = 20, n_dimred = 20)
    
    set.seed(100)
    sce.sel = runUMAP(sce.sel, use_dimred="MNN", perplexity = 20, n_dimred = 20)
    
    
    plotTSNE(sce.sel, colour_by="cluster", size_by = "total_features_by_counts") + fontsize + ggtitle("seurat - graph base clustering")
    
    plotUMAP(sce.sel, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
      fontsize + ggtitle("seurat -- graph based clustering")
    
    
    
  }
  
  dev.off()
  
}

Find.Gene.Markers.with.scran = FALSE
if(Find.Gene.Markers.with.scran){
  markers <- findMarkers(sce, my.clusters)
  
  marker.set <- markers[["1"]]
  head(marker.set, 10)
  
  top.markers <- rownames(marker.set)[marker.set$Top <= 5]
  plotHeatmap(sce, features=top.markers, 
              columns=order(sce$cluster), 
              colour_columns_by=c("cluster"),
              cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5)) 
  
  
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







##################################################
##################################################
## Project: transcriptional priming mechanism by tbx in C.elegans
## Script purpose: analysis the single-cell RNA-seq data
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Mon Feb 19 14:43:38 2018
##################################################
##################################################
version.DATA = 'scRNA_R6875_R7116_R7130_R7130redo_R7133_R7926'
version.analysis =  paste0(version.DATA, '_20190703')

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
Merge.technicalRep = TRUE

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
if(length(grep("out_gene.featureCounts.txt", colnames(aa)))>0) {
  aa = aa[, -grep("out_gene.featureCounts.txt", colnames(aa))]
}

##########################################
# manually make design matrix 
##########################################
if(Manually.Specify.sampleInfos.filtering.4scRNAseq){
  library(stringi)
  library("openxlsx")
  
  design = data.frame(samples = colnames(aa)[-1], stringsAsFactors = FALSE)
  design = data.frame(sapply(design, function(x) gsub('[.]', '_', x)), stringsAsFactors = FALSE)
  
  design$flowcell.lane = sapply(design$samples, function(x) paste0(unlist(strsplit(as.character(x), "_"))[c(1:2)], collapse = "_"))
  design$sampleIDs = sapply(design$samples, function(x) unlist(strsplit(as.character(x), "_"))[3])
  #design$sampleIDs = gsub('_\\w+$', '', design$sampleIDs)
  design$barcodes = sapply(design$samples, function(x) unlist(strsplit(as.character(x), "_"))[4])
  
  ##########################################
  # here manually 
  ##########################################
  design$request = NA;
  design$request[which(design$flowcell.lane == "CCVTBANXX_8")] = "R6875"
  design$request[which(design$flowcell.lane == "CCVBPANXX_1")] = "R7116"
  design$request[which(design$flowcell.lane == "HHG5KBGX9_1")] = "R7130"
  design$request[which(design$flowcell.lane == "HHGHNBGX9_1")] = "R7130"
  
  design$request[which(design$flowcell.lane == "HLWTCBGX9_1")] = "R7130"
  design$request[which(design$flowcell.lane == "CCYTEANXX_4")] = "R7130"
  design$request[which(design$flowcell.lane == "CD2GTANXX_5")] = "R7133"
  design$request[which(design$flowcell.lane == "H7KNYBGXB_1")] = "R7926"
  

  design$seqInfos = paste0(design$request, "_", design$flowcell.lane)
  #jj = grep("CCVTBANXX_8.76090_", colnames(aa))
  #aa = aa[, c(1, jj)]
  #colnames(aa)[-1] = sapply(colnames(aa)[-1], function(x) gsub("CCVTBANXX_8.76090_", "", x))
  kk = which(!is.na(design$request))
  design = design[kk, ]
  aa = aa[, c(1, kk+1)]
  
}

##########################################
# aggregated quality controls from nf-RNAseq 
##########################################
if(Aggregate.nf.QCs.plots.in.designMatrix){
  #load(file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_RNA_seq.Rdata'))
  
  source("scRNAseq_functions.R")
  dirs.all = c('../../../Ariane/R7116_R7130_scrnaseq/results_all/multiqc_data_1', 
               "../../../Ariane/R6875_scRNAseq/results_all_3rd/MultiQC/multiqc_data", 
               "../../R7130_redo_R7133/results_v2/multiqc_data_1", 
               "../../R7926/results_v2/multiqc_data")
  QCs.nf = aggrate.nf.QCs(dirs.all)
  QCs.nf$Sample = gsub("#", "_", QCs.nf$Sample)
  
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
library(SingleCellExperiment)
library(scater)
library(scRNA.seq.funcs)
library(scran)
options(stringsAsFactors = FALSE)

Manual.vs.outlier.filtering = FALSE
load(file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_design_sampleInfos_QCs_nf_RNA_seq.Rdata'))

source("scRNAseq_functions.R")
counts = convertGeneNames.forCountTable(aa)
gg.Mt = find.particular.geneSet("Mt")
gg.ribo = find.particular.geneSet("ribo")

##########################################
# compare tehnical replicates,
# merge them 
# or benchmark batch correction methods
##########################################
if(Merge.technicalRep){
  
  source("scRNAseq_functions.R")
  pdfname = paste0(resDir, "/scRNAseq_check_technicalRep.pdf")
  pdf(pdfname, width=12, height = 6)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  xx = merge.techinical.replicates(design = design, counts = counts, 
                                   sampleInfos.techRep = list(c("R7130_HHG5KBGX9_1", "R7130_HLWTCBGX9_1"), 
                                                              c("R7130_HHGHNBGX9_1", "R7130_CCYTEANXX_4", "R7133_CD2GTANXX_5")))
  dev.off()
  
  xx1 = xx$design
  xx2 = xx$counts
  
  design = xx1
  counts = xx2
  
  save(design, counts, file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_design_technicalRepMerged.Rdata'))

}

##########################################
# Import SingleCellExperiment and scater packages for the QC and table cleaning
# several steps will be proceded:
# 1) general overview of data quality: sequencing depth, mapping rate, assignment rate, rRAN codamination for each sequencing lane
# 2) clean the cells 
# 3) clean genes
##########################################
source("scRNAseq_functions.R")
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)

load(file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_design_technicalRepMerged.Rdata'))

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
pdf(pdfname, width=12, height = 6)
par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

# some general statistics for each request and lane
plotColData(sce, y="pct_counts_Ribo", x="seqInfos") + ggtitle("% of rRNA contamination")
plotColData(sce, y="pct_counts_Mt", x="seqInfos") + ggtitle("% of Mt")

plotColData(sce, y="log10_total_counts", x="seqInfos") + ggtitle("total nb of reads mapped to transcripts")

plotColData(sce, y="total_features_by_counts", x="seqInfos") + ggtitle("total nb of genes")

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
filter_by_MT = sce$pct_counts_Mt < 5
table(filter_by_MT)

sce$use <- (
    filter_by_expr_features & # sufficient features (genes)
    filter_by_total_counts & # sufficient molecules counted
    # filter_by_ERCC & # sufficient endogenous RNA
    filter_by_MT # remove cells with unusual number of reads in MT genes
)
table(sce$use)


## comparison between manual filtering and automatic ouliter filtering
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

####################
## filter lowly expressed (and probably too highly expressed genes)
####################
load(file=paste0(RdataDir, version.DATA, '_QCed_cells_filtered_SCE.Rdata'))

pdfname = paste0(resDir, "/scRNAseq_QCs_genes_filterting_", version.analysis, ".pdf")
pdf(pdfname, width=16, height = 10)
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
# Section : scRNA-seq data normalization 
# codes original from Hemberg's course
# https://hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#data-visualization
# Among all normalization methods, DESeq2 and scran normalization will be used
# consider seurat for UMI counts, but here we are just working on the read counts
# Many normalization have been proposed
# Here we test main two methods: TMM from edgeR or from DESeq2
# and the one from scran
# the PCA and some other plots were used to assess the normalization
# After testing different normalization using technical replicates, cpm, DESeq2 and scran yield similar results and scran is slightly better 
# than DESeq2 which is better than cpm 
# !!! parameters of scran normalization required extra care; double check the size factor from scran vs library size
# the current parameters for scran: 
# qclust <- quickCluster(sce.qc, min.size = 100, method = 'igraph')
# sce.qc <- computeSumFactors(sce.qc, clusters = qclust)
########################################################
########################################################
load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_SCE.Rdata'))
library(scRNA.seq.funcs)
library(scater)
library(scran)
options(stringsAsFactors = FALSE)

Normalization.Testing = FALSE

reducedDim(sce) <- NULL
endog_genes <- !rowData(sce)$is_feature_control

if(Normalization.Testing){
  source("scRNAseq_functions.R")
  
  pdfname = paste0(resDir, "/scRNAseq_filtered_normalization_testing.pdf")
  pdf(pdfname, width=14, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  test.normalization(sce, Methods.Normalization = c("cpm", "DESeq2", "scran"), min.size = 100)
  
  dev.off()
}

##########################################
# select normalization method: scran
##########################################
set.seed(1000)    
clusters <- quickCluster(sce, min.size = 100, method="igraph")
table(clusters)

sce <- computeSumFactors(sce, clusters = clusters)

## quick check for size factors calculated by scran
summary(sizeFactors(sce))
plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
     xlab="Library size (millions)", ylab="Size factor", pch=16)


sce <- normalize(sce, exprs_values = "counts", return_log = TRUE)

save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))




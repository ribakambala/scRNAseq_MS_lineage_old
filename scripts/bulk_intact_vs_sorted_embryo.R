##########################################################################
##########################################################################
# Project: Aleks and Ariane MS lineage with scRNAseq project
# Script purpose: to compare the intact and sorted c elegans embryo 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Nov  5 12:18:12 2018
##########################################################################
##########################################################################
library("openxlsx")
require('DESeq2')
source('scRNAseq_functions.R')

### data verision and analysis version   
version.Data = 'bulk_sc_rnaseq_Rxxxx_R6875_all_v1';
version.analysis = paste0("_", version.Data, "_20181105")

### Directories to save results
design.file = "../exp_design/Bulkseq_Rxxxx.xlsx"
data.file = "../data/R6875_Rxxx_merged_gene_counts.txt"

resDir = "../results/bulk_intact_vs_sorted"
tabDir =  paste0(resDir, "/tables/")
RdataDir = paste0(resDir, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

Processing.design.matrix = TRUE
##################################################
##################################################
## Section: import design matrix and prepare count table
##################################################
##################################################
design = read.xlsx(design.file, sheet = 1, colNames = TRUE)
#design = design[which(design$stage != "L1.early"),]

if(Processing.design.matrix){
  xx = data.frame(design, stringsAsFactors = FALSE);
  design = data.frame(xx$Sample.ID, xx$Brief.Sample.Description, xx$Additional.Info, xx$Genetic.Background, stringsAsFactors = FALSE)
  colnames(design) = c('SampleID',  'treatment', 'adaptor.concentration', 'genotype')
  
  kk = grep("dissociated", design$treatment)
  design$treatment[kk] = "dissociated.sorted"

  kk = c(1:6)
  design$adaptor.concentration[kk] = "Tn5.1.75"
  design$adaptor.concentration[-kk] = "Tn5.1.150"
  design$genotype = "N2"
  #design = design[order(design$stage, design$treatment), ]
  #design$tissue.cell[which(design$genotype=="henn-1_mutant" & design$promoter=="no_promoter")] = "whole.body_no_promoter"
}

## make the data table
xlist = data.file
#xlist<-list.files(path=paste0(dataDir), pattern = "*.txt", full.names = TRUE) ## list of data set to merge
#spikes.file = xlist[grep("spikeIns_", xlist)]
#spikes.file = spikes.file[grep("_old", spikes.file, invert = TRUE)]
#xlist = xlist[grep("cel_", xlist)]

if(length(xlist)>1){
  all = NULL
  ggs = NULL
  #counts = list()
  for(n in 1:length(xlist)){
    ff = read.delim(xlist[n], sep='\t', header = TRUE);
    colnames(ff)[1] = "gene";
    
    if(n == 1){
      all = data.frame(ff, stringsAsFactors = FALSE);
      colnames(all)[1] = "gene"
    }else{
      ggs = ff[,1];
      gg.union = unique(union(all$gene, ggs))
      all = data.frame(gg.union, 
                       all[match(gg.union, all$gene), -1],
                       ff[match(gg.union, ggs), -1], 
                       stringsAsFactors = FALSE)
      colnames(all)[1] = "gene"
    }
  }
}else{
  all = read.delim(xlist, sep = "\t", header = TRUE)
}


# processing count table
all = process.countTable(all=all, design = design);

#### Import Sample information and table of read counts
#design = read.xlsx(paste0(design.file), sheet = 1, colNames = TRUE)
#colnames(design)[c(1:2)] = c("SampleID", "genotype")
#design$genotype[grep("WT", design$genotype)] = "WT"
#design$tissue.cell[which(design$genotype=="henn-1_mutant" & design$promoter=="no_promoter")] = "whole_animal_no_promoter"

#all = read.delim(data.file, sep = "\t", header = TRUE)
# processing count table
#all = process.countTable(all=all, design = design);

spikes = read.delim(spikes.file, sep="\t", header = TRUE, row.names = 1)
spikes = t(as.matrix(spikes))
spikes = data.frame(gene=rownames(spikes), spikes, stringsAsFactors = FALSE)
spikes = process.countTable(all=spikes, design = design, select.Total.count = FALSE)

all = rbind(spikes, all);
#spikes = process.countTable(spikes, design)

save(design, all, file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))

##################################################
##################################################
## Section: spike-in normalization and QC for cpm and spike-in normalization
##################################################
##################################################
load(file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))
read.count = all[, -1];
kk = c(1:nrow(design))

design.matrix = data.frame(sample=colnames(read.count)[kk], design[kk, ])
raw = as.matrix(read.count[,kk])
raw[which(is.na(raw))] = 0
### start enrichment analysis 
raw = floor(raw)
rownames(raw) = all$gene

####################
## QC for cpm 
####################
QC.for.cpm = FALSE
if(QC.for.cpm){
  #treat = length(unique(design$treatment[kk]));
  #index.qc = c(3, 5)[which(c(length(unique(design.matrix$genotype)), length(unique(design.matrix$promoter)))>1)]
  index.qc = c(1, 3,4)
  
  source("RNAseq_Quality_Controls.R")
  pdfname = paste0(resDir, "/Data_qulity_assessment_early_L1_and_all", version.analysis, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  Check.RNAseq.Quality(read.count=read.count[, kk], design.matrix = design.matrix[, index.qc])
  dev.off()
}

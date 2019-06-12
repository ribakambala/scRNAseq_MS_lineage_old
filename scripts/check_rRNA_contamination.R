##########################################################################
##########################################################################
# Project: Aleks and Ariane MS-lineage 
# Script purpose: check the potential reson for rRNA contamination
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Nov  5 17:09:33 2018
##########################################################################
##########################################################################
stat.file = "../data/R6875_Rxxx_multiqc_general_stats.txt"
plate.file = "../exp_design/R6875_position_info_in_plate_pool_tags.csv"

stat = read.delim(stat.file, sep = "\t", header = TRUE, as.is = c(1))
positions = read.csv(plate.file, header = TRUE)

kk = grep("76090", stat$Sample)
stat = stat[kk, ]

resDir = "../results/rRNA_contamination"
if(!dir.exists(resDir)){dir.create(resDir)}

pdfname = paste0(resDir, "/check_potential_reasons_for_rRNA_contamination.pdf")
pdf(pdfname, width = 12, height = 10)

par(mfrow=c(1, 1))
hist(stat$Biotype.Counts_mqc.generalstats.biotype_counts.percent_rRNA, breaks = 12, xlab = "% of rRNA")

source("RNAseq_Quality_Controls.R")
yy = stat[, -1]

par(mfrow=c(3, 4))
for(n in 1:ncol(yy)) {
  plot(yy[, c(11, n)], xlab=colnames(yy)[11], ylab = colnames(yy)[n])
}

#pairs(yy, lower.panel=NULL, upper.panel=panel.fitting, main = "pairwise-comparisons")

stat$nb_row = NA;
stat$nb_col = NA;
adaptor_bcs = paste0(positions$adaptor_tag, positions$adaptor_secondary_tag)
for(n in 1:nrow(stat))
{
  #n = 1
  bc = gsub("CCVTBANXX_8#76090_", "", stat$Sample[n])
  kk = which(adaptor_bcs==bc)
  stat$nb_col[n] = positions$adaptor_number[kk]
  stat$nb_row[n] = positions$adaptor_secondary_number[kk]
}

index.col = unique(stat$nb_col)
index.col = index.col[order(index.col)]
index.row = unique(stat$nb_row)
index.row = index.row[order(index.row)]

yy = matrix(NA, nrow = length(unique(stat$nb_row)), ncol=length(unique(stat$nb_col)))

for(n in 1:nrow(stat)){
  #n = 1
  yy[which(index.row==stat$nb_row[n]), which(index.col==stat$nb_col[n])] = stat$Biotype.Counts_mqc.generalstats.biotype_counts.percent_rRNA[n]
}

colnames(yy) = index.col
rownames(yy) = index.row

colors <- colorRampPalette( (brewer.pal(9, "Blues")) )(255)
pheatmap(yy, cluster_rows = FALSE, cluster_cols = FALSE,
         col = colors)


dev.off()
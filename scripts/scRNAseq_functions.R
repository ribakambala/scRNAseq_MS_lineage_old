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

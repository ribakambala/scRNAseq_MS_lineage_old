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


## common functions for normalizations from Hemberg lab
calc_cpm <- function (expr_mat, spikes = NULL) 
{
  norm_factor <- colSums(expr_mat[-spikes, ])
  return(t(t(expr_mat)/norm_factor)) * 10^6
}
calc_sf <- function (expr_mat, spikes = NULL) 
{
  geomeans <- exp(rowMeans(log(expr_mat[-spikes, ])))
  SF <- function(cnts) {
    median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 
                              0)])
  }
  norm_factor <- apply(expr_mat[-spikes, ], 2, SF)
  return(t(t(expr_mat)/norm_factor))
}
calc_uq <- function (expr_mat, spikes = NULL) 
{
  UQ <- function(x) {
    quantile(x[x > 0], 0.75)
  }
  uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
  norm_factor <- uq/median(uq)
  return(t(t(expr_mat)/norm_factor))
}
calc_cell_RLE <-  function (expr_mat, spikes = NULL) 
{
  RLE_gene <- function(x) {
    if (median(unlist(x)) > 0) {
      log((x + 1)/(median(unlist(x)) + 1))/log(2)
    }
    else {
      rep(NA, times = length(x))
    }
  }
  if (!is.null(spikes)) {
    RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
  }
  else {
    RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
  }
  cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
  return(cell_RLE)
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

##########################################################################
##########################################################################
# Project: single-cell RNAseq data
# Script purpose: aggregate the quality controls from nf-RNAseq 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jan 15 15:18:33 2019
##########################################################################
##########################################################################
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
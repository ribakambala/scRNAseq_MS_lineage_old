##########################################################################
##########################################################################
# Project: single-cell RNAseq data
# Script purpose: aggregate the quality controls from nf-RNAseq 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jan 15 15:18:33 2019
##########################################################################
##########################################################################
aggrate.nf.QCs = function(dir)
{
  # dir = '../../../Ariane/R7116_R7130_scrnaseq/results_all/'
  folds = paste0(dir, c("STAR", "rseqc", "featureCounts"))
  
  for (f in folds)
  {
    f = folds[1]
    if(dir.exists(f)){
      print(f);
      if(basename(f) == "STAR"){
        files = list.files(path = paste0(f, "/logs"), pattern = "*.final.out", full.names = TRUE)
        
          
        for(n in 1:length(files))
        {
          if(n %%100 ==0) cat('---', n,  '---\n')
          xx = read.delim(files[n], header = FALSE, as.is = 1)
          xx = xx[which(xx[, 1]) == "Number of input reads" |
                    ]
        }
        
      } 
        
      
      
    }
  }
}
##########################################################################
##########################################################################
# Project: Aleks' single cell RNA-seq 
# Script purpose: iterative clustering function
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct 16 11:38:07 2019
##########################################################################
##########################################################################
cal_cor_sd_time = function(var, nb.timepoints = 9){
  # var = xx[4, ]
  nb.vec = length(var) / nb.timepoints
  dv = t(matrix(var, nrow = nb.vec, byrow = TRUE))
  dv.cor = cor(dv)
  dv.sd = apply(dv, 2, sd)
  
  return(c(cor.min = min(dv.cor, na.rm = TRUE),  cor.median = median(dv.cor, na.rm = TRUE),
           sd.min = min(dv.sd, na.rm = TRUE), sd.median = median(dv.sd, na.rm = TRUE)))
    
}

cal_autocor_stationaryTest_sd = function(var, plot = FALSE)
{
  # var = bulk[19, ]; plot = TRUE
  ac = acf(var, lag.max = 7, plot = plot)
  if(plot){ plot(var) }
  ac = abs(ac$acf[-1])
  #acf(as.numeric(bulk[1, ]), lag.max = 10, plot = TRUE)
  ac.test = adf.test(var) # p-value < 0.05 indicates the TS is stationary
  # kpss.test(tsData)var
  sds = sd(var,  na.rm = TRUE)
  
  return(c(ac.max =  max(ac), pval.ac = ac.test$p.value, sd = sds))
  
}


Test.timer.genes.with.HashimshonyLineageData = function(timers)
{
  
}

find.timer.genes = function(Test.timer.genes.with.HashimshonyData = FALSE)
{
   dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
   require('openxlsx')
   
   double.check.lineage.independent.genes = FALSE
   if(double.check.lineage.independent.genes){
     lineages = read.xlsx(paste0(dataDir.Hashimsholy, "/Hashimshony_datasets_lineage_11timepoints.xlsx"), sheet = 1, startRow = 9, colNames = FALSE)
     
     names = lineages[c(1, 2), ]
     xx = lineages[-c(1, 2), ]
     
     lineages.names = c('AB_T', 'MS_T', 'E_T', "C_T", 'P3_T')
     timepoints = c(1:11)
     colnames(xx) = c('ensEMBL.IDs', 'Wormbase',  as.vector(t(outer(lineages.names, timepoints, paste, sep="")))) 
     
     geneMapping = xx[, c(1, 2)] 
     
     #save(geneMapping, file = paste0(dataDir.Hashimsholy, "/geneNames_mapping_Wormbase.Rdata"))
     
     rownames(xx) = xx[, 1]
     xx = xx[, -c(1, 2)]
     kk = grep('T1$|T2', colnames(xx))
     xx = xx[, -kk]
     xx = as.matrix(xx)
     ss = apply(xx, 1, sum)
     jj = which(ss>2^5)
     
     xx = xx[jj, ]
     
     xx = log2(xx + 2^-8)
     
     yy = t(apply(xx, 1, cal_cor_sd_time))
     
     plot(yy[, 1], yy[, 3], cex = 0.5)
     plot(yy[, 2], yy[, 4], cex = 0.5)
   }else{
    load( file = paste0(dataDir.Hashimsholy, "/geneNames_mapping_Wormbase.Rdata")) 
   }
   
   bulk = read.delim(file = paste0(dataDir.Hashimsholy, "/GSE50548_Whole_embryo_interval_timecourse.txt"), 
                    row.names = 1,  header = TRUE, sep = "\t")
   
   mm = match(rownames(bulk), geneMapping$ensEMBL.IDs)
   rownames(bulk) = geneMapping$Wormbase[mm]
   
   bulk = as.matrix(bulk)
   ss = apply(bulk, 1, sum)
   kk = which(ss>2^5)
   bulk = bulk[kk, ]
   
   bulk = log2(bulk + 2^-6)
   
   library(tseries)
   
   yy = t(apply(bulk, 1, cal_autocor_stationaryTest_sd))
   
   plot(yy[, c(1, 3)], cex = 0.5);
   abline(v = 0.6, col='red')
   abline(h = 1.5, col = 'red')
   
   yy = data.frame(yy)
   sels = which(yy$ac.max>0.6 & yy$sd > 1.5)
   
   timers = bulk[sels, ]
   
   plot(timers[100,], type='b')
   
   if(Test.timer.genes.with.HashimshonyData){
     Test.timer.genes.with.HashimshonyLineageData(timers)
   }
   
   return(timers)
}





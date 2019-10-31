##########################################################################
##########################################################################
# Project: Aleks' single cell RNA-seq 
# Script purpose: iterative clustering function
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct 16 11:38:07 2019
##########################################################################
##########################################################################
cal_cor_sd_time = function(vec, nb.timepoints = 9){
  # vec = xx[4, ]
  nb.vec = length(vec) / nb.timepoints
  dv = t(matrix(vec, nrow = nb.vec, byrow = TRUE))
  dv.cor = cor(dv)
  dv.sd = apply(dv, 2, sd)
  
  return(c(cor.min = min(dv.cor, na.rm = TRUE),  cor.median = median(dv.cor, na.rm = TRUE),
           sd.min = min(dv.sd, na.rm = TRUE), sd.median = median(dv.sd, na.rm = TRUE)))
    
}

cal_autocor_stationaryTest_sd = function(vec, plot = FALSE, lag.length = 15)
{
  library(tseries)
  options(warn=-1)
  # vec = bulk[2369, ]; plot = TRUE
  ac = acf(vec, lag.max = lag.length, plot = plot)
  ac = abs(ac$acf[-1])
  #acf(as.numeric(bulk[1, ]), lag.max = 10, plot = TRUE)
  adf.test = adf.test(vec) # p-value < 0.05 indicates the TS is stationary
  box.test = Box.test(vec, lag=lag.length, type="Ljung-Box")
  # kpss.test(tsData)vec
  sds = sd(vec,  na.rm = TRUE)
  
  if(plot){ plot(vec, main = paste0('sd = ', sds)) }
  
  return(c(ac.max =  max(ac), sd = sds, pval.box = box.test$p.value, pval.adf = adf.test$p.value))
  
}

Test.timingEstimate.with.HashimshonyLineages = function()
{
  require('openxlsx')
  lineages = read.xlsx(paste0(dataDir.Hashimsholy, "/Hashimshony_datasets_lineage_11timepoints.xlsx"), sheet = 1, startRow = 9, colNames = FALSE)
  
  names = lineages[c(1, 2), ]
  xx = lineages[-c(1, 2), ]
  
  lineages.names = c('AB_T', 'MS_T', 'E_T', "C_T", 'P3_T')
  timepoints = c(1:11)
  colnames(xx) = c('ensEMBL.IDs', 'Wormbase',  as.vector(t(outer(lineages.names, timepoints, paste, sep="")))) 
  
  ## keep genes considered as timers
  xx = xx[match(rownames(timers), xx[,2]), ]
  rownames(xx) = xx[, 2]
  
  
  ## import another table from GEO
  aa = read.table(paste0(dataDir.Hashimsholy, "/GSE50548_Blastomere_developmental_event_timecourse.txt"), header = TRUE)
  aa = aa[match(xx$ensEMBL.IDs, rownames(aa)), ]
  rownames(aa) = rownames(xx)
  test = mapply(aa, FUN=as.numeric)
  rownames(test) = rownames(aa)
  test = log2(test + 2^-6)
  
  xx = xx[, -c(1, 2)]
  test = mapply(xx, FUN=as.numeric)
  rownames(test) = rownames(xx)
  test = log2(test + 2^-6)
  
  estimate.timing(test[, 1], timers)
  
  
  for(kk in c(1:55)){
    # kk = 55
    
  }
  
  
  
}

find.timer.genes = function(plot.test = FALSE)
{
  dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
  load(file = paste0(dataDir.Hashimsholy, "/annotMapping_ensID_Wormbase_GeneName.Rdata"))
  # annot = read.csv(file = "../../../../annotations/BioMart_WBcel235_noFilters.csv", header = TRUE)
  #geneMapping$Gene.name = annot$Gene.name[match(geneMapping$Wormbase, annot$Gene.stable.ID)]
  #save(geneMapping, file = paste0(dataDir.Hashimsholy, "/annotMapping_ensID_Wormbase_GeneName.Rdata"))
  
  bulk = read.delim(file = paste0(dataDir.Hashimsholy, "/GSE50548_Whole_embryo_interval_timecourse.txt"), 
                    row.names = 1,  header = TRUE, sep = "\t")
  
  bulk = as.matrix(bulk)
  ss = apply(bulk, 1, sum)
  
  ## filter the genes by
  kk = which(ss>0)
  bulk = bulk[kk, ]
  ss = ss[kk]
  
  kk = which(ss>10^2)
  bulk = bulk[kk, ]
  
  bulk = log2(bulk + 2^-6)
  
  ## the main result: max of autocorrection, sd, pval for Ljung-Box test for independence test (pval < 0.05 non-stationary time series), 
  ## Augmented Dickey–Fuller (ADF) t-statistic test for unit root test(pval large for non-stationary time sereis) 
  ## one reference web: https://rpubs.com/richkt/269797
  timers = t(apply(bulk, 1, cal_autocor_stationaryTest_sd))
  timers = data.frame(timers)
  
  if(plot.test){
    
    par(mfrow = c(1,2))
    plot(-log10(timers$pval.box), timers$ac.max, cex = 0.2)
    abline(h= 0.5, col='red', lwd= 2.0)
    abline(v = 3, col ='red', lwd = 2.0)
    plot(-log10(timers$pval.box), timers$sd, cex = 0.2)
    abline(v = 3, col ='red', lwd = 2.0)
    abline(h= 1, col='red', lwd= 2.0)
    
    #plot(timers[100,], type='b')
  }
  
  mm = match(rownames(timers), geneMapping$ensEMBL.IDs)
  geneNames = geneMapping$Gene.name[mm]
  jj = !is.na(geneNames)
  timers = timers[jj, ]
  rownames(timers) = geneNames[jj]
  
  # save(timers, file =paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval.Rdata"))
  
  return(timers)
  
}


estimate.timing.with.timer.genes = function(vec, timerGenes.pval=0.001, timerGenes.ac=0.5, timerGenes.sd = 1.5, loess.span = 1.5, reprocess.timer.genes = FALSE)
{
  # input: 
  # vec -- a vector of gene expression with names of WBGeneID 
  
  vec = test[, kk]
  corrs = rep(NA, length = ncol(timers))
  names(corrs) = colnames(timers)
  timepoints = gsub('mpfc_','',  colnames(timers))
  timepoints[c(1:8)] = c('0', '5', '10', '20', '40', '50', '60', '70')
  timepoints = as.numeric(timepoints)
  
  for(n in 1:ncol(timers)) {
    #ii = which(vec > (-6) & timers[,n] > (-6))
    ii = which(timers[,n] > (-6))
    # ii = c(1:nrow(timers))
    corrs[n] = cor(vec[ii], timers[ii,n], method = 'pearson')
  }
  predFun = loess(corrs ~ timepoints, span = 0.15)
  
  prediction = predict(predFun, timepoints, se = FALSE)
  index = which(prediction==max(prediction))
  estimate.timing = timepoints[index]
  
  plot(timepoints, corrs)
  points(timepoints, prediction, type = 'l', lwd=1.5, col='darkblue')
  points(timepoints[index], prediction[index], col = 'darkred', cex = 2.0, pch = 16)
  
  cat(colnames(test)[kk], " vs. ", estimate.timing, "min -- corr : ", prediction[index], "\n")
  
  return(estimate.timing)
  
}









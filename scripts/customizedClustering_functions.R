##########################################################################
##########################################################################
# Project: Aleks' single cell RNA-seq 
# Script purpose: iterative clustering function
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct 16 11:38:07 2019
##########################################################################
##########################################################################
mycor = function(x, y){
  x = as.numeric(x)
  y = as.numeric(y)
  if(sd(x) < 10^(-6) | sd(y) <10^(-6)){
    return(NA)
  }else{
    return(cor(x, y, method = 'pearson'))
  }
}

cal_cor_sd_time = function(vec, ref, nb.timepoints = 8, Plot.test = FALSE){
  # vec = xx[10, ]; ref = refs[10, ];nb.timepoints = 8
  nb.vec = length(vec) / nb.timepoints
  dv = matrix(vec, nrow = nb.vec, byrow = TRUE)
 
  dv.cor = apply(dv, 1, function(x) mycor(x, ref))
  #dv.sd = apply(dv, 2, sd)
  names(dv.cor) = c('AB', "MS", "E", 'C')
  
  if(Plot.test){
    matplot(c(1:nb.timepoints), t(dv))
    points(c(1:nb.timepoints), ref, type ='l', lwd =2.0, col='darkblue')
    
  }
  
  return(dv.cor)
    
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

##########################################
# this function is to refine timerGenes by adding the correlation between AB, MS, E, C and whole-embryo
# to further select genes independent lineages
# even though the evidence are limited but this is the best we can do
##########################################
refine_timerGenes_using_lineages = function(timers)
{
  ## import the timer genes
  Modify.timeponts.Remove.1cellStage.forTimers = FALSE
  if(Modify.timeponts.Remove.1cellStage.forTimers){
    dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
    load(file = paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval_plus_timepoints.Rdata"))
    
    kk.1cell = which(colnames(timers) == 'X1_cell')
    if(length(kk.1cell) == 1) 
    {
      timers = timers[, -kk.1cell]
      timepoints = timepoints[-1]
      timepoints = timepoints + 20
      timepoints[1] = 0
    }
    save(timers, timepoints, file = paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval_plus_timepoints.Rdata"))
  }else{
    dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
    load(file = paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval_plus_timepoints.Rdata"))
  }
  
  ## import the table of five lineages from Hashimsholy et al. paper
  ## import the table downloaded from GEO
  dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
  load(file = paste0(dataDir.Hashimsholy, "/annotMapping_ensID_Wormbase_GeneName.Rdata"))
  
  aa = read.table(paste0(dataDir.Hashimsholy, "/GSE50548_Blastomere_developmental_event_timecourse.txt"), header = TRUE)
 
  Compare.lineage.tables.from.Wormbase.vs.GEO = FALSE
  if(Compare.lineage.tables.from.Wormbase.vs.GEO){
    require('openxlsx')
    lineages = read.xlsx(paste0(dataDir.Hashimsholy, "/Hashimshony_datasets_lineage_11timepoints.xlsx"), sheet = 1, startRow = 9, colNames = FALSE)
    
    mm = match(rownames(aa), lineages$X1)
    lineages = lineages[mm, ]
    
    #names = lineages[c(1, 2), ]
    aa2 = lineages[, -c(1, 2)]
    rownames(aa2) = lineages$X1
    lineages.names = c('AB_T', 'MS_T', 'E_T', "C_T", 'P3_T')
    colnames(aa2) = c(as.vector(t(outer(lineages.names, c(1:11), paste, sep=""))))
    
    ## keep genes considered as timers
    # xx = xx[match(rownames(timers), xx[,2]), ]
    # rownames(xx) = xx[, 2]
    # xx = xx[, -c(1, 2)]
    xx = mapply(aa2, FUN=as.numeric)
    rownames(xx) = rownames(aa2)
    xx = log2(xx + 2^-6)
    
    test = mapply(aa, FUN=as.numeric)
    rownames(test) = rownames(aa)
    test = log2(test + 2^-6)
    
    n = 4
    plot(test[, n], xx[, n], cex = 0.4);
    abline(0, 1, lwd=2.0, col ='red')
  }
  
  names = geneMapping$Gene.name[match(rownames(aa), geneMapping$ensEMBL.IDs)]
  jj = !is.na(names)
  bb = aa[jj, ]
  rownames(bb) = names[jj]
  
  test = mapply(bb, FUN=as.numeric)
  rownames(test) = rownames(bb)
  test = log2(test + 2^-6)
  
  test = test[match(rownames(timers), rownames(test)), ]
  ## select the corresponding time points for timer genes and lineages 
  test = test[, c(1:44)]
  #kk = grep('tr|E8...|E8..', colnames(test))
  xx = test[, -c(9:11, (9:11) +11, (9:11) +22, (9:11) +33)]
  
  refs = timers[, -c(1:4)]
  time.sels = c(1, 2, 3, 4, 7,8,9,12)
  cat(timepoints[time.sels], "\n")
  
  refs = refs[, time.sels]
  
  temporal.cor = matrix(NA, ncol = 4, nrow = nrow(refs))
  for(n in 1:nrow(refs)){
    cat(n, "\n")
    temporal.cor[n, ] = cal_cor_sd_time(xx[n, ], refs[n, ], nb.timepoints = 8)
  }
 
  colnames(temporal.cor) = c('AB', "MS", "E", 'C')
  tcors = data.frame(temporal.cor)
  
  rownames(tcors) = rownames(timers)
  
  save(timers, timepoints, tcors, file = paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval_plus_timepoints_tempCorrelation.Rdata"))
  
  cutoffs = 0.6
  kk = which(tcors$AB>cutoffs & tcors$MS>cutoffs & tcors$E>cutoffs & tcors$C>cutoffs & timers$pval.box<0.0001)
  length(kk)
  
  pdfname = paste0("../results/clustering_combining_variousInfos/timerGenes_tempCorrelation_lineages_vs_wholeEmbry.pdf")
  pdf(pdfname, width=10, height = 6)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(n in 1:100){
    cal_cor_sd_time(xx[kk[n], ], refs[kk[n], ], nb.timepoints = 8, Plot.test = TRUE)
  }
  
  dev.off()
  #plot(-log10(timers$pval.box), tcors$AB, cex =0.5)
  
}

##########################################
# the main function to find genes showing clear time-depenent patterns in the whole embry data
##########################################
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
  
  timers = data.frame(timers, bulk)
  
  mm = match(rownames(timers), geneMapping$ensEMBL.IDs)
  geneNames = geneMapping$Gene.name[mm]
  jj = !is.na(geneNames)
  timers = timers[jj, ]
  rownames(timers) = geneNames[jj]
  
  # save(timers, file =paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval.Rdata"))
  
  return(timers)
}


estimate.timing.with.timer.genes = function(vec, timerGenes.pval=0.001, timerGenes.ac=0.5, timerGenes.sd = 0, loess.span = 0.15,
                                            use = 'lowFilter.timers', lowFilter.threshold.timer = -6, lowFilter.threshold.target = -6, 
                                            reprocess.timer.genes = FALSE, PLOT.test = FALSE)
{
  # input: 
  # vec -- a vector of gene expression with names of WBGeneID 
  
  # vec = test[, kk]; reprocess.timer.genes = FALSE;  use = 'lowFilter.timers'; loess.span = 0.15;
  if(reprocess.timer.genes){
    find.timer.genes(plot.test = FALSE)
  }else{
    dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
    load(file =paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval.Rdata"))
    
  }
  
  ## select the timer genes using pval and autocorrelation
  sels.timerGenes = which(timers$ac.max > timerGenes.ac & timers$pval.box < timerGenes.pval)
  timers = timers[sels.timerGenes, -c(1:4)]
  
  # select the common genes shared by timers and vec
  sharedGenes = intersect(rownames(timers), names(vec))
  vec = vec[match(sharedGenes, names(vec))]
  timers = timers[match(sharedGenes, rownames(timers)),]
  
  #vec = test[, kk]
  corrs = rep(NA, length = ncol(timers))
  names(corrs) = colnames(timers)
  timepoints = gsub('mpfc_','',  colnames(timers))
  timepoints[c(1:8)] = c('0', '5', '10', '20', '40', '50', '60', '70')
  timepoints = as.numeric(timepoints)
  
  for(n in 1:ncol(timers)) {
    if(use == 'lowFilter.timers') ii = which(timers[,n] > lowFilter.threshold.timer)
    if(use ==  'lowFilter.target') ii = which(vec > lowFilter.threshold.target) 
    if(use ==  'lowFilter.both') ii = which(vec > lowFilter.threshold.target & timers[,n] > lowFilter.threshold.timer)
    if(use ==  'all') ii = c(1:nrow(timers))
    #cat('using timer genes -- ', length(ii), "\n")
    corrs[n] = cor(vec[ii], timers[ii,n], method = 'pearson')
  }
  
  predFun = loess(corrs ~ timepoints, span = loess.span)
  prediction = predict(predFun, timepoints, se = FALSE)
  estimate.timing = timepoints[which(prediction==max(prediction))[1]]
  
  if(PLOT.test){
    par(mfrow = c(1,1))
    index = which(prediction==max(prediction))
    plot(timepoints, corrs)
    points(timepoints, prediction, type = 'l', lwd=1.5, col='darkblue')
    points(timepoints[index], prediction[index], col = 'darkred', cex = 2.0, pch = 16)
  }
  
  return(estimate.timing)
  
}

## here is a fast version of estimate.timing.with.timer.genes function
## where use = 'lowFilter.target'
fast.estimate.timing.with.timer.genes = function(vec, timers, timepoints, timerGenes.pval=0.001, timerGenes.ac=0.5, timerGenes.sd = 0, loess.span = 0.15,
                                            use = 'lowFilter.target', lowFilter.threshold.target = 0, 
                                            reprocess.timer.genes = FALSE, PLOT.test = FALSE)
{
  if(reprocess.timer.genes){
    timers = find.timer.genes(plot.test = FALSE)
    
    timepoints = gsub('mpfc_','',  colnames(timers)[-c(1:4)])
    timepoints[c(1:8)] = c('0', '5', '10', '20', '40', '50', '60', '70')
    timepoints = as.numeric(timepoints)
    
    save(timers, timepoints, file = paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval_plus_timepoints.Rdata"))
  }
  
  ## select the timer genes using pval and autocorrelation
  sel.vec = which(vec > lowFilter.threshold.target) 
  vec = vec[sel.vec]
  
  sels.timerGenes = which(timers$ac.max > timerGenes.ac & timers$pval.box < timerGenes.pval)
  timers = timers[sels.timerGenes, -c(1:4)]
  
  # select the common genes shared by timers and vec
  sharedGenes = intersect(rownames(timers), names(vec))
  # cat('nb of timer genes used :', length(sharedGenes), "\n")
  vec = vec[match(sharedGenes, names(vec))]
  timers = timers[match(sharedGenes, rownames(timers)),]
  
  corrs = apply(timers, 2, function(x) cor(vec, x, method = 'pearson'))
  
  predFun = loess(corrs ~ timepoints, span = loess.span)
  prediction = predict(predFun, timepoints, se = FALSE)
  estimate.timing = timepoints[which(prediction==max(prediction))[1]]
  
  if(PLOT.test){
    par(mfrow = c(1,1))
    index = which(prediction==max(prediction))
    plot(timepoints, corrs)
    points(timepoints, prediction, type = 'l', lwd=1.5, col='darkblue')
    points(timepoints[index], prediction[index], col = 'darkred', cex = 2.0, pch = 16)
  }
  
  return(estimate.timing)
  
}

### test timing.estimation using lineages from Hashimshony et al. paper
Test.timingEstimate.with.HashimshonyLineages = function(timerGenes.pval = 0.001, loess.span = 0.15, use = 'lowFilter.target')
{
  ## import the table downloaded from GEO
  dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
  aa = read.table(paste0(dataDir.Hashimsholy, "/GSE50548_Blastomere_developmental_event_timecourse.txt"), header = TRUE)
  load(file = paste0(dataDir.Hashimsholy, "/annotMapping_ensID_Wormbase_GeneName.Rdata"))
  
  names = geneMapping$Gene.name[match(rownames(aa), geneMapping$ensEMBL.IDs)]
  jj = !is.na(names)
  bb = aa[jj, ]
  rownames(bb) = names[jj]
  
  test = mapply(bb, FUN=as.numeric)
  rownames(test) = rownames(bb)
  test = log2(test + 2^-6)
  
  estimation = c()
  for(kk in c(1:55)){
    # kk = 55
    timing = estimate.timing.with.timer.genes(vec = test[,kk], PLOT.test = TRUE, timerGenes.pval=timerGenes.pval, use = use, 
                                              loess.span = loess.span)
    cat(colnames(test)[kk], " vs. ", timing, "min \n")
    estimation = c(estimation, timing)
  }
  
  experiments = rep(c(0, 20, 40, 60, 90, 110, 140, 180, 200, 300, 400), 5)
  
  plot(experiments, estimation, type='n', main = paste0('timerGene.pval = ', signif(timerGenes.pval, d= 3), ' & loess.span = ', signif(loess.span, d=2)))
  for(n in c(1:5)) {
    index = seq(11*(n-1)+1, 11*n, by = 1)
    points(experiments[index], estimation[index], pch = n, cex = 1.2, col = 'darkblue')
  }
  abline(0, 1, lwd=2.0, col='darkred')
  abline(-20, 1, lwd=1.0, col='darkred', lty = 2)
  abline(-40, 1, lwd=1.0, col='darkred', lty= 2)
  abline(-60, 1, lwd=1.0, col='darkred', lty= 2)
  legend('topleft', legend = c('AB', 'MS', 'E', 'C', 'P3'), pch = c(1:5), col = 'darkblue', bty = 'n')
  
}


sc.estimateTiming.with.timer.genes = function(sce, fastEstimate = TRUE, timerGenes.pval = 0.001, loess.span = 0.15, lowFilter.threshold.target = 5, 
                                              PLOT.test = FALSE)
{
  ## extract the gene expression matrix from sce object
  test = logcounts(sce)
  
  if(fastEstimate){
    
    # timerGenes.pval = 0.001; loess.span = 0.15; lowFilter.threshold.target = 5
    dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
    load(file = paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval_plus_timepoints.Rdata"))
    estimation.fast = rep(NA, ncol(test))
    for(kk in c(1:ncol(test))){
      estimation.fast[kk] = fast.estimate.timing.with.timer.genes(vec = test[,kk], timers = timers, timepoints = timepoints,
                                                                  PLOT.test = PLOT.test,
                                                                  timerGenes.pval= timerGenes.pval, loess.span = loess.span, 
                                                                  lowFilter.threshold.target = lowFilter.threshold.target)
    }
    sce$timing.est = estimation.fast
    
  }else{
    cat('For the sake of speed, use fastEstimate \n')
    # estimation = apply(test, 2, estimate.timing.with.timer.genes, PLOT.test = PLOT.test, 
    #                    timerGenes.pval= timerGenes.pval, loess.span = loess.span, 
    #                    use = use, lowFilter.threshold.target = lowFilter.threshold.target)
  }
  
  return(sce)
}

Test.sc.estimateTiming.with.timer.genes = function(sce)
{
  
  ## extract the gene expression matrix from sce object
  test = logcounts(sce)
  set.seed(2019)
  index.test = sample(c(1:ncol(test)), 500)
  #index.test = c(1:ncol(test))
  test = test[, index.test]
  
  if(fastEstimate){
    
    # timerGenes.pval = 0.0001; loess.span = 0.2; lowFilter.threshold.target = 5;  PLOT.test = FALSE
    dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
    load(file = paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval_plus_timepoints.Rdata"))
    
    
    library(tictoc)
    
    tic('for loop ')
    estimation.fast = rep(NA, ncol(test))
    estimation.test = rep(NA, ncol(test))
    for(kk in c(1:ncol(test))){
      # kk = 132; PLOT.test = TRUE; cat(estimation.fast[kk], "-- ", estimation.test[kk], "\n")
      if(kk%%200 == 0) cat(kk, "\n")
      estimation.fast[kk] = fast.estimate.timing.with.timer.genes(vec = test[,kk], timers = timers, timepoints = timepoints,
                                                                  PLOT.test = PLOT.test,
                                                                  timerGenes.pval= timerGenes.pval, loess.span = loess.span, 
                                                                  lowFilter.threshold.target = lowFilter.threshold.target)
      
      estimation.test[kk] = fast.estimate.timing.with.timer.genes(vec = test[,kk], timers = timers, timepoints = timepoints,
                                                                  PLOT.test = PLOT.test,
                                                                  timerGenes.pval= timerGenes.pval, loess.span = loess.span, 
                                                                  lowFilter.threshold.target = 6)
      
    }
    toc()
    
    diffs = abs(estimation.fast - estimation.test)
    
    sce$timing.est = estimation.fast
    
    Test = FALSE
    if(Test){
      par(mfrow = c(1, 3))
      plot(sce$FSC_log2[index.test], estimation.fast, type='p', 
           main = paste0('timerGene.pval = ', signif(timerGenes.pval, d= 3), ' & loess.span = ', signif(loess.span, d=2)))
      plot(sce$BSC_log2[index.test], estimation.fast, type='p') 
      plot(sce$FSC_log2[index.test], sce$BSC_log2[index.test], type = 'p', cex = 0.5)
      
    }
    
  }else{
    tic()
    estimation = rep(NA, ncol(test))
    for(kk in c(1:ncol(test))){
      #cat(kk, "\n")
      estimation[kk] = estimate.timing.with.timer.genes(vec = test[,kk], PLOT.test = FALSE, 
                                                        timerGenes.pval= timerGenes.pval, loess.span = loess.span, 
                                                        use = use, lowFilter.threshold.target = 2)
      
    }
    toc()
    
    tic()
    estimation = apply(test, 2, estimate.timing.with.timer.genes, PLOT.test = FALSE, 
                       timerGenes.pval= timerGenes.pval, loess.span = loess.span, 
                       use = use, lowFilter.threshold.target = 2)
    toc()
  }
  
  
  sce.test0 = sc.estimateTiming.with.timer.genes(sce, fastEstimate = TRUE, timerGenes.pval = 0.0001, loess.span = 0.5, lowFilter.threshold.target = 5, 
                                                 PLOT.test = FALSE)
  sce.test7 = sc.estimateTiming.with.timer.genes(sce, fastEstimate = TRUE, timerGenes.pval = 0.0001, loess.span = 0.15, lowFilter.threshold.target = 5, 
                                                 PLOT.test = FALSE)
  
  par(mfrow = c(1, 1))
  plot(sce.test0$timing.est, sce.test7$timing.est); 
  abline(0,1,col = 'red'); 
  abline(-30, 1, col = 'red', lty =2);  
  abline(30, 1, col = 'red', lty = 2)
  abline(-60, 1, col = 'red', lty =3);  
  abline(60, 1, col = 'red', lty = 3)
  
  par(mfrow = c(1, 3))
  plot(sce$FSC_log2, sce.test0$timing.est, type='p', cex = 0.5)
  plot(sce$BSC_log2, sce.test0$timing.est, type='p', cex = 0.5) 
  plot(sce$FSC_log2, sce$BSC_log2, type = 'p', cex = 0.5)
  
}


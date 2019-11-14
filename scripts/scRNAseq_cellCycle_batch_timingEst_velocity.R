########################################################
########################################################
# Section : Dealing with batch confouner
# first thing to try is the MNN method from scran pacakge after contacting with the first 
# author, Laleh Haghverdi and Aaron Lun
# !!!!!! to correct the batch, the main ideal is to 
# 1) select highly variable genes 2) using these genes to reduce the dimensions and correct the batch in the low dimension
# 3) the corrected output will be low dimensional output, which will be used to define the distance for clustering and also trajectory
# or the corrected gene expression matrix (but not used for DE analysis), from which PCs can also be found by PCA
# Taken together, HGV, dimensionality reduction, batch correction are together. 
# cluster and DE analysis will be the next two steps. 
# in addition, the DE analysis will be performed with original counts by providing cluster labels and batches(
# http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/de.html#2_blocking_on_uninteresting_factors_of_variation)
########################################################
########################################################
path2AleksFolder = '/Volumes/groups/cochella/Aleks/bioinformatics/GitHub/scRNAseq_MS_lineage'
version.DATA = 'scRNA_8613_full'
version.analysis =  paste0(version.DATA, '_20191029')

dataDir = paste0("../data/")
resDir = paste0("../results/", version.analysis)
tabDir = paste0("../results/", version.analysis, "/tables/")
RdataDir = paste0("../results/", version.analysis, "/Rdata/")
RdataDirfromAleks = paste0(path2AleksFolder, "/results/", version.analysis, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

##########################################
# Remove the cell cycle confounder 
# here we choose to use Seurat to regress out the cell cycle effect
# we need to train the cells to identify the cell cycle phase
# this could be more complicated than expected
##########################################
correct.cellCycle = FALSE
Import.processed.data.from.Aleks = TRUE
if(Import.processed.data.from.Aleks){
  load(file=paste0(RdataDirfromAleks, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))
}else{
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))
}

if(correct.cellCycle){
  source("scRNAseq_functions.R")
  # cellCycle.correction(sce, method = "seurat")
  #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected.Rdata'))
  #sce_old = sce
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2.Rdata')) 
  
  library(scater)
  p1 = scater::plotPCA(
    sce,
    run_args = list(exprs_values = "logcounts", feature_set = c(1:500)), 
    size_by = "total_counts",
    #size_by = "total_features_by_counts",
    colour_by = "seqInfos"
  ) + ggtitle(paste0("PCA -- "))
  
  p2 = scater::plotPCA(
    sce,
    run_args = list(exprs_values = "logcounts_seurat", feature_set = c(1:500)), 
    size_by = "total_counts",
    #size_by = "total_features_by_counts",
    colour_by = "seqInfos"
  ) + ggtitle(paste0("PCA -- "))
  
  multiplot(p1, p2, cols = 2)
}else{
  save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2.Rdata'))
}

##########################################
##########################################
# Prepare the data for batch correction (fastMNN from scran used here) and clustering 
# Feature selection (HVGs) is performed
# here we just use the scran trandVar()
# further test to follow for other methods
# (https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby011/4898116)
# But many of them, e.g. BASics, Brennecke works better with spike-in 
##########################################
##########################################

#load(file = file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos.Rdata'))
facsInfo.Added = TRUE

##########################################
# 1) add the facs information in the metadata
# 2) estimate timing for celes using timer genes
# 3) have better time estimation by combing both information
##########################################
if(Import.processed.data.from.Aleks){
  load(file = paste0(RdataDirfromAleks, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrectedv2_facsInfos.Rdata'))
  save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrectedv2_facsInfos.Rdata')) 
}else{
  if(facsInfo.Added){
    load(file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos.Rdata'))
  }else{
    load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2.Rdata')) 
    
    source("scRNAseq_functions.R")
    sce = Integrate.FACS.Information(sce)
    
    save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrectedv2_facsInfos.Rdata')) 
  }
}

# logcounts(sce) =  assay(sce, "logcounts_seurat_SG2MCorrected")
table(sce$nb.cells)
sce = sce[, which(sce$nb.cells == 1)]

sce$FSC_log2 = 3/2*log2(sce$FSC)
sce$BSC_log2 = 3/2*log2(sce$BSC)

plotColData(sce,
            x = "FSC_log2",
            y = "BSC_log2",
            colour_by = "request",
            point_size = 1
)

##########################################
# here estimat the timing with timer genes
### first test 5 lineages from Hashimsholy et al. paper
##########################################
reEstimate.timing.using.timer.genes = FALSE
if(reEstimate.timing.using.timer.genes){
  Test.Hashimshony_lineages = FALSE
  
  if(Test.Hashimshony_lineages){
    pdfname = paste0("../results/clustering_combining_variousInfos/test_timing_estimation_Hashimshony_lineages_test_with_improvedTimerGenes_v2.pdf")
    pdf(pdfname, width=10, height = 6)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    source('customizedClustering_functions.R')
    Test.timingEstimate.with.HashimshonyLineages(fastEstimate = TRUE, timerGenes.pval = 0.0001, lineageCorrs = 0.7,  loess.span = 0.5, lowFilter.threshold.target = 5, 
                                                 PLOT.test = FALSE)
    
    dev.off()
    
  }
  
  ## Here we are sampling a range of parameters and timing estimation were done with each of them
  ## Whereby we assess the sensibility of our timingEst
  ## this will take some time to finish
  source('customizedClustering_functions.R')
  
  timingEst = c()
  for(pv in c(0.001, 0.0001, 0.00001))
  {
    for(cutoff.expr in c(4, 5, 6))
    {
      for(s in c(0.3, 0.5, 0.7))
      {
        cat('pv = ', pv, ' cutoff.expr = ', cutoff.expr, 's = ', s, "\n")
        sce.test = sc.estimateTiming.with.timer.genes(sce, fastEstimate = TRUE, timerGenes.pval = pv, lineageCorrs = 0.5, loess.span = s, 
                                                      lowFilter.threshold.target = cutoff.expr)
        timingEst = rbind(timingEst, sce.test$timing.est)
      }
    }
  }
  
  #save(sce, timingEst, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrectedv2_facsInfos_timeEst_tmp.Rdata')) 
  
  timingEst = t(timingEst)
  timingEst = as.matrix(timingEst)
  #colnames(timingEst) = c(as.vector(t(outer(c(0.001, 0.0001, 0.00001), c(4:6), paste, sep=""))))
  find.one.close.to.mean = function(x){
    # x = timingEst[1, ]
    difs = abs(x - mean(x))
    return(x[which(difs == min(difs))[1]])
  }
  sce$timingEst = apply(timingEst, 1, find.one.close.to.mean)
  sce$timingEst.sd = apply(timingEst, 1, sd)
  
  save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos_timingEst.Rdata')) 
}else{
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos_timingEst.Rdata'))
}


par(mfrow = c(1, 3))
plot(sce$FSC_log2, sce$timingEst, type='p', cex = 0.5)
plot(sce$BSC_log2, sce$timingEst, type='p', cex = 0.5) 
plot(sce$FSC_log2, sce$BSC_log2, type = 'p', cex = 0.5)

cdata = colData(sce)
cdata = data.frame(cdata[, c((ncol(cdata)-3): ncol(cdata))])
#cdata$timingEst = cdata$timingEst/50

cdata$timing.group = NA
cdata$timing.group[which(cdata$timingEst < 50)] = 1
cdata$timing.group[which(cdata$timingEst >= 450)] = 10
for(n in 2:9){cdata$timing.group[which(cdata$timingEst >= (n-1)*50 & cdata$timingEst < n*50)] = n}

cdata$timing.sd.group = 3
cdata$timing.sd.group[which(cdata$timingEst.sd<30)] = 1
cdata$timing.sd.group[which(cdata$timingEst.sd>=30 & cdata$timingEst.sd<60)] = 2
cdata$timing.sd.group = as.factor(cdata$timing.sd.group)

ggplot(cdata, aes(x=FSC_log2, y=BSC_log2, color=timingEst, shape = timing.sd.group)) +
  geom_point() + 
  scale_color_gradientn(colours = rainbow(10))

sce$timingEst = as.factor(sce$timingEst)
sce$timingEst.group = as.factor(cdata$timing.group)
sce$timingEst.sd.group = as.factor(cdata$timing.sd.group)

save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos_timingEstGroups.Rdata'))
plotColData(sce,
            x = "FSC_log2",
            y = "BSC_log2",
            colour_by = "timingEst.group",
            point_size = 1
)
########################################################
########################################################
# Section : Batch Correction
# Batch correction using fastMNN from scran
# here we are using fastMNN
########################################################
########################################################

##########################################
# Feature selection (select HGVs) for batch corrections
# there are two options: batch-specific or use batch as block
# https://www.bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html
##########################################
library(scRNA.seq.funcs)
library(scater)
library(SingleCellExperiment)
library(scran)
library(kBET)
set.seed(1234567)
options(stringsAsFactors = FALSE)

load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos_timingEstGroups.Rdata'))

sce$timingEst = as.factor(sce$timingEst)
sce$timingEst.group = as.factor(sce$timingEst.group)

# choose the batches (either plates or request)
# here we choose the request as batch
batches = sce$request 
bc.uniq = unique(batches)
sce$batches <- batches

Use.fastMNN = TRUE
Norm.Vars.per.batch = TRUE # HVGs for each batch or not 
Rescale.Batches = FALSE # scaling data in each batch or not 
k.mnn = 20
cos.norm = TRUE
nb.pcs = 50

batch.sequence.to.merge = c('R7130', 'R8612', 'R8526', 'R7926', # 3 plates for each request
                            'R6875','R7116','R8613','R8348') # 1 plate for each request

order2correct = match(batch.sequence.to.merge, bc.uniq) 
  #c(12, 13, # the same 
  #                10, 9, 8, 5, 6, 7,  8, 11, 1, 2, 3, 4)
#order2correct = c(15,14,13,12,11,10,9,8, 1, 2, 3, 4, 5,6,7 )
#order2correct = c(3, 4, 1, 2)

## double chekc  the mering order in the batch correction
source('customizedClustering_functions.R')
kk = match(sce$request, c('R7130', 'R7926', 'R8526', 'R8612'))
kk = match(sce$request, c('R7130', 'R6875', 'R7116', 'R8348'))
plotColData(sce[,which(!is.na(kk))],
            x = "FSC_log2",
            y = "BSC_log2",
            colour_by = "request",
            point_size = 1
)

#cat("merging order for batch correction :\n", paste0(bc.uniq[order2correct], collapse = "\n"), "\n")
for(n in 1:length(order2correct)){
  
  #n = 11
  kk = order2correct[n]
  
  p = plotColData(sce[, which(sce$batches== bc.uniq[kk])],
              x = "FSC_log2",
              y = "BSC_log2",
              colour_by = "timingEst",
              point_size = 1
              
  )
  plot(p)
  
  cat('#', n, 'batch:',  bc.uniq[kk], ': ', length(which(sce$batches == bc.uniq[kk])), 'cells \n')
  
}


source("scRNAseq_functions.R")
HVGs = find.HVGs(sce, Norm.Vars.per.batch = Norm.Vars.per.batch, method = "scran", ntop = 2000)
gene.chosen = match(HVGs, rownames(sce))

cat("nb of HGV : ", length(gene.chosen), "\n")


if(Use.fastMNN){
  ## rescaling for each batch is recommended by the author
  ## We adjust the size factors with multiBatchNorm() to make them more comparable across batches. 
  ## This mitigates differences in scale and variance in the log-expression values between batches, especially between technologies.
  ## https://www.bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/multibatch.html#3_feature_selection_across_batches
  if(Rescale.Batches){
    ttxt = c("nout = multiBatchNorm(")
    for(n in 1:length(bc.uniq)){
      if(n != length(bc.uniq)) {
        ttxt = paste0(ttxt, "sce[, which(sce$batches == '", bc.uniq[n], "')], ")
      } else {
        ttxt = paste0(ttxt, "sce[, which(sce$batches == '", bc.uniq[n], "')])")
      }
    }
    eval(parse(text = ttxt)) ## rescaling done
    
    kk2check = 1
    par(mfrow=c(1,1))
    plot(sizeFactors(nout[[kk2check]]), sizeFactors(sce[, which(sce$batches == bc.uniq[kk2check])])); abline(0, 1, col='red')
  }
  
  original = list()
  fscs = c()
  #original0 = list()
  for(n in 1:length(bc.uniq)){
    #xx = nout[[n]];
    if(Rescale.Batches){
      original[[n]] = logcounts((nout[[n]][gene.chosen, ]))
    }else{
      original[[n]] = logcounts((sce[gene.chosen, which(sce$batches == bc.uniq[n])])) 
    }
    fscs = c(fscs, sce$FSC_log2[which(sce$batches == bc.uniq[n])])
  }
  
  # Slightly convoluted call to avoid re-writing code later.
  # Equivalent to fastMNN(GSE81076, GSE85241, k=20, d=50, approximate=TRUE)
  # Comments from Aaron: https://www.bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html#5_examining_the_effect_of_correction
  # The k= parameter specifies the number of nearest neighbours to consider when defining MNN pairs. This should be interpreted as the minimum frequency of each cell type or state in each batch.
  # Larger values will improve the precision of the correction by increasing the number of MNN pairs. It also provides some robustness to violations of the assumption that the batch vector is orthogonal to the biological subspace (Haghverdi et al. 2018), by allowing the neighbour search to ignore biological variation in each batch to identify the correct MNN pairs.
  # However, larger values of k can also reduce accuracy by allowing incorrect MNN pairs to form between cells of different types. Thus, we suggest starting with the default k and increasing it if one is confident that the same cell types are not adequately merged across batches. This is better than starting with a large k as incorrect merging is much harder to diagnose than insufficient merging.
  # When BSPARAM=IrlbaParam(deferred=TRUE), fastMNN() uses methods from the irlba package to perform the principal components analysis quickly. While the run-to-run differences should be minimal, it does mean that set.seed() is required to obtain fully reproducible results. The deferred= argument instructs fastMNN() to sacrifice some numerical precision for greater speed.
  set.seed(1001)
  mnn.out <- do.call(fastMNN, c(original, list(k=k.mnn, cos.norm = cos.norm, d=nb.pcs, auto.order=order2correct,
                                               approximate=TRUE)))
  dim(mnn.out$corrected)
  mnn.out$batch
  Rle(mnn.out$batch)
  #metadata(mnn.out)$merge.info$pairs[[1]]
  reducedDim(sce, "MNN") <- mnn.out$corrected
  sce$mnn_Batch <- as.character(mnn.out$batch)
  sce
  
  ##########################################
  # check the effectiveness of batch correction with MNN
  # 1) check the MNN pairs and lost variances
  # 2) visualization wiht PCA, and UMAP before and after correction
  # 3) kBET test
  ##########################################
  pdfname = paste0(resDir, "/scRNAseq_filtered_test_MNNbatchCorrection_effectiveness_v1.pdf")
  pdf(pdfname, width=14, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  set.seed(1001)
  nb.MNN.to.use = 20
  sce <- runUMAP(sce, use_dimred="MNN", n_dimred = nb.MNN.to.use, ncomponents = 2, scale_features = FALSE,
                 method = c("umap-learn"))
  
  p = plotUMAP(sce, ncomponents = 2, colour_by="timingEst.group", size_by = "FSC_log2", point_size= 0.01) + ggtitle("Corrected")
  plot(p)
  
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
  p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  multiplot(p1, p2, cols = 2)
  
  sce = runTSNE(sce, use_dimred="MNN", perplexity = 30, n_dimred = nb.MNN.to.use)
  p = plotTSNE(sce, colour_by="timingEst", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  
  p1 = plotTSNE(sce, colour_by="pha-4", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  p2 = plotTSNE(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  multiplot(p1, p2, cols = 2)
  
  # mnn.out$pairs
  source("scRNAseq_functions.R")
  Check.MNN.pairs(mnn.out, fscs)
  
  #omat <- do.call(cbind, original)
  #sce.qc <- SingleCellExperiment(list(logcounts=omat))
  set.seed(1000)
  with.var <- do.call(fastMNN, c(original,
                                 list(k=k.mnn, cos.norm = cos.norm, d=nb.pcs, auto.order=order2correct, approximate=TRUE,
                                      compute.variances=TRUE)))
  with.var$lost.var
  
  fsc.boundary = c(27.2, 28.6)
  bsc.boundary = c(24.5, 26.8)
  plotColData(sce,
              x = "FSC_log2",
              y = "BSC_log2",
              colour_by = "mnn_Batch",
              shape_by = "mnn_Batch"
              
  ) + geom_vline(xintercept = fsc.boundary, linetype="dashed", color = "blue", size=0.5) +
    geom_hline(yintercept= bsc.boundary, linetype="dashed", color = "red", size=0.5) 
  
  sel.tmp = which(sce$FSC_log2 > fsc.boundary[1] & sce$FSC_log2 < fsc.boundary[2] 
                  & sce$BSC_log2 > bsc.boundary[1] & sce$BSC_log2 < bsc.boundary[2])
  #sce.tmp = sce[gene.chosen, which(sce$mnn_Batch > 2)]
  sce.tmp = sce[, sel.tmp]
  #sce.tmp = sce
  dim(sce.tmp)
  
  sce.tmp <- runPCA(sce.tmp, ncomponents = 50, ntop=nrow(sce.tmp), method="irlba", exprs_values = "logcounts", scale_features = TRUE)
  
  plotPCA(sce.tmp, colour_by="batches") + ggtitle("Original")
  
  dff = as.data.frame(reducedDim(sce.tmp, "MNN"))
  colnames(dff) = paste0("PC", c(1:ncol(dff)))
  dff$batches = sce.tmp$batches
  ggp = ggplot(data=dff, aes(PC1, PC2, color=batches)) + geom_point(size=2) 
  plot(ggp);
  
  # Using irlba to set up the t-SNE, for speed.
  #set.seed(100)
  #osce <- runTSNE(sce.tmp, use_dimred="PCA", perplexity = 20)
  #ot <- plotTSNE(osce, colour_by="mnn_Batch") + ggtitle("Original")
  #set.seed(100)
  #csce <- runTSNE(sce, use_dimred="MNN", perplexity = 20)
  #ct <- plotTSNE(csce, colour_by="mnn_Batch") + ggtitle("Corrected")
  #multiplot(ot, ct, cols=2)
  
  set.seed(100)
  osce <- runUMAP(sce.tmp, use_dimred="PCA", n_dimred = nb.MNN.to.use)
  ot <- plotUMAP(osce, colour_by="mnn_Batch", size_by = "FSC_log2") + ggtitle("Original")
  set.seed(100)
  csce <- runUMAP(sce.tmp, use_dimred="MNN", n_dimred = nb.MNN.to.use)
  ct <- plotUMAP(csce, colour_by="mnn_Batch", size_by = "FSC_log2") + ggtitle("Corrected")
  multiplot(ot, ct, cols=2)
  
  set.seed(100)
  osce <- runTSNE(sce.tmp, use_dimred="PCA", perplexity = 20, n_dimred = nb.MNN.to.use)
  ot <- plotTSNE(osce, colour_by="mnn_Batch", size_by = "FSC_log2") + ggtitle("Original")
  set.seed(100)
  csce <- runTSNE(sce.tmp, use_dimred="MNN", perplexity = 20, n_dimred = nb.MNN.to.use)
  ct <- plotTSNE(csce, colour_by="mnn_Batch", size_by = "FSC_log2") + ggtitle("Corrected")
  multiplot(ot, ct, cols=2)
  
  kbet.orig <- kBET(
    df = t(reducedDim(sce.tmp, "PCA")), 
    batch = sce.tmp$mnn_Batch,
    heuristic = TRUE,
    do.pca = FALSE,
    verbose = TRUE, 
    addTest = FALSE,
    n_repeat = 200,
    plot = TRUE)
  
  kbet.bc = kBET(
    df = t(reducedDim(sce.tmp, "MNN")), 
    batch = sce.tmp$mnn_Batch,
    do.pca = FALSE,
    heuristic = FALSE,
    verbose = TRUE, 
    addTest = FALSE,
    n_repeat = 200,
    plot = TRUE)
  
  dev.off()
  
}

save(sce, file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_bcMNN.Rdata'))

##########################################
# Convert sce object to Seurat object 
# check UMAP and tSNE 
##########################################
require(Seurat)
pbmc = as.Seurat(sce)

#pbmc = Seurat::RunPCA(pbmc, pbmc, )
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#pbmc <- RunPCA(pbmc, features = HVGs)

pbmc <- Seurat::RunUMAP(pbmc, dims = 1:15, reduction = "MNN", 
                        reduction.key = "umap", n.neighbors = 20, repulsion.strength = 1)

DimPlot(pbmc, reduction = "umap", group.by = 'timingEst')

pbmc <- Seurat::RunUMAP(pbmc, dims = 1:15, reduction = "pca", 
                        reduction.key = "umap.pca", n.neighbors = 20, repulsion.strength = 1)

DimPlot(pbmc, reduction = "umap", group.by = 'timingEst')


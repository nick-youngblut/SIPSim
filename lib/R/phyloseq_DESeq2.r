#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: phyloseq_DESeq2.r [options] <phyloseq>

options:
  <phyloseq>   Phyloseq object file.
  --log2=<l>   Log2 fold change cutoff.
               [Default: 0.25]
  --hypo=<h>   altHypothesis tested by DESeq
               ("greaterAbs","greater","less")
               [Default: NULL]
  --cont=<c>   Control libraries. (comma-separated list)
               [Default: 1]
  --treat=<t>  Treatment libraries. (comma-separated list)
               [Default: 2]
  -h           Help' -> doc

opts = docopt(doc)
cont = unlist(strsplit(opts[['--cont']], split=','))
treat = unlist(strsplit(opts[['--treat']], split=','))

# packages
pkgs <- c('phyloseq', 'DESeq2')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# main
log2.cut = as.numeric(opts[['--log2']])
## import
if(opts[['<phyloseq>']] == '-'){
  con = pipe("cat", "rb")
  physeq = readRDS(con)
} else {
  physeq = readRDS(opts[['<phyloseq>']])
}

## library must be a character
physeq.sd = sample_data(physeq)
physeq.sd$library = as.factor(as.character(physeq.sd$library))
physeq.sd$condition = physeq.sd$library
sample_data(physeq) = physeq.sd

## making deseq2 object
dds = phyloseq_to_deseq2(physeq, ~condition)

### making 'condition' factor
condition = as.character(dds$condition)
for (x in cont){
  condition = ifelse(condition == x, 'control', condition)
}
for (x in treat){
  condition = ifelse(condition == x, 'treatment', condition)
}
dds$condition = factor(condition, levels=c('control', 'treatment'))

if(any(is.na(dds$condition))){
  stop('\nSome control/treatments are NA. Check --cont & --treat.\n')
}


# calculate geometric means prior to estimate size factors
## This method is not sensitive to zeros
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans=geoMeans)

## DESeq
dds = DESeq(dds, fitType='local', test='Wald')
if(opts[['--hypo']] != 'NULL'){
  res = results(dds,
    lfcThreshold=log2.cut,
    altHypothesis=opts[['--hypo']],
    independentFiltering=FALSE)  
} else {
  res = results(dds, independentFiltering=FALSE)
}


## p-adjust as in Pepe-Ranney et al. (2015)
beta = res$log2FoldChange
betaSE = res$lfcSE
res$p = pnorm(beta, log2.cut, betaSE, lower.tail = FALSE)
res$padj.BH = p.adjust(res$p, "BH")

## writing
con = pipe("cat", "wb")
saveRDS(res, con)

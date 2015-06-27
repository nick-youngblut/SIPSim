#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: phyloseq_DESeq2.r [options] <phyloseq>

options:
  <phyloseq>   Phyloseq object file.
  -h           Help' -> doc

opts = docopt(doc)


# packages
pkgs <- c('phyloseq', 'DESeq2')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# main
## import
if(opts[['<phyloseq>']] == '-'){
  con = pipe("cat", "rb")
  physeq = readRDS(con)
} else {
  physeq = readRDS(opts[['<phyloseq>']])
}

## library must be a character
physeq.sd = sample_data(physeq)
physeq.sd$library = as.character(physeq.sd$library)
sample_data(physeq) = physeq.sd

## making deseq2 object
dds = phyloseq_to_deseq2(physeq, ~library)
dds$library = relevel(as.factor(dds$library), 1)

# calculate geometric means prior to estimate size factors
## This method is not sensitive to zeros
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans=geoMeans)

## DESeq
dds = DESeq(dds, fitType='local')
res = results(dds)

## writing
con = pipe("cat", "wb")
saveRDS(res, con)

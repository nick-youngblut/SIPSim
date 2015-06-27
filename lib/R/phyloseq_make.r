#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: phyloseq_make.r [options] <OTU>

options:
  <OTU>     OTU table in "wide" format.
  -s=<s>    Sample table 
  -h        Help' -> doc

opts = docopt(doc)


# packages
pkgs <- c('phyloseq')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}


# main
## OTU
otu = read.delim(opts[['<OTU>']], check.names=F, row.names=1)
otu = otu_table(otu, taxa_are_rows=T)

## Sample
if(! is.null(opts[['-s']])){
  samp = read.delim(opts[['-s']], row.names=1)
  samp = sample_data(samp)
} else {
  samp = FALSE
}

## phyloseq
physeq = phyloseq(otu, samp)

## writing
con = pipe("cat", "wb")
saveRDS(physeq, con)

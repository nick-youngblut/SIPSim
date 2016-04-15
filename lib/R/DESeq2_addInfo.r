#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: DESeq2_rare-dominant.r [options] <DESeq2> <comm> <BD_shift>

options:
  <DESeq2>         DESeq2 results object file ("-" if from STDIN).
  <comm>           Community file.
  <BD_shift>       BD-shift file.
  -h               Help
description:
  Adding information to DESeq object file:
    * taxon abundance
    * amount of BD distribution shift
  
  Taxon abundance is defined as the mean percent
  abundance over all communities in the community file.

  Output: DESeq table written to STDOUT
' -> doc

opts = docopt(doc)

# packages
pkgs <- c('dplyr')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# options
abund.cut = as.numeric(opts[['-c']])

# main
## import
if(opts[['<DESeq2>']] == '-'){
  con = pipe("cat", "rb")
#  deseq.res = readRDS(con)
  deseq.res = read.delim(con, sep='\t')
} else {
 # deseq.res = readRDS(opts[['<DESeq2>']])
  deseq.res = read.delim(opts[['<DESeq2>']], sep='\t')
}
x = suppressPackageStartupMessages(as.data.frame(deseq.res))

comm = read.delim(opts[['comm']], sep='\t')
BD.shift = read.delim(opts[['BD_shift']], sep='\t')


## edit
### DESeq2 object
deseq.res = as.data.frame(deseq.res)
deseq.res$taxon = rownames(deseq.res)
### Community file
comm = comm %>%
  group_by(taxon_name) %>%
    summarize(mean_rel_abund_perc = mean(rel_abund_perc)) %>%
      ungroup() %>%
        mutate(taxon_name = as.character(taxon_name),
               rank = dense_rank(mean_rel_abund_perc))
### BDshift
BD.shift = BD.shift %>%
      filter(lib2 == 2) %>%
        mutate(taxon = as.character(taxon))
### table join
tbl.j = inner_join(deseq.res, BD.shift, c('taxon' = 'taxon'))
tbl.j = inner_join(tbl.j, comm, c('taxon' = 'taxon_name'))
rownames(tbl.j) = as.character(tbl.j$taxon)

## parsing and writing
## writing
con = pipe("cat", "wb")
saveRDS(tbl.j, con)

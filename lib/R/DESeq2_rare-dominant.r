#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: DESeq2_rare-dominant.r [options] <DESeq2> <comm>

options:
  <DESeq2>         DESeq2 results object file ("-" if from STDIN).
  <comm>           Community file used to parse DESeq2 results
                   into "dominant" and "rare" taxa.
  -o=<o>           Basename for output files.
                   [Default: DESeq2]
  -c=<c>           Abundance cutoff (percent abundance) for defining
                   rare and dominant taxa.
                   [Default: 1]
  -h               Help
description:
  Parse DESeq2 output into abundance-class bins.

  This can be used to see if DESeq2 accuracy is
  biased by taxon abundance.

  Taxon abundance is defined as the mean percent
  abundance over all communities in the community file.
' -> doc

opts = docopt(doc)
padj.cut = as.numeric(opts[['--padj']])

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
#  deseq.res = readRDS(opts[['<DESeq2>']])
  deseq.res = read.delim(opts[['<DESeq2>']], sep='\t')  
}
x = suppressPackageStartupMessages(as.data.frame(deseq.res))
comm = read.delim(opts[['comm']], sep='\t')


## edit
### DESeq2 object
deseq.res = as.data.frame(deseq.res)
deseq.res$taxon = rownames(deseq.res)
### Community file
comm = comm %>%
  group_by(taxon_name) %>%
    summarize(mean_rel_abund_perc = mean(rel_abund_perc)) %>%
      ungroup() %>%
        mutate(taxon_name = as.character(taxon_name))
### table join
tbl.j = inner_join(deseq.res, comm, c('taxon' = 'taxon_name'))
rownames(tbl.j) = tbl.j$taxon

## parsing and writing
### rare
tbl.j.rare = tbl.j %>% filter(mean_rel_abund_perc <= abund.cut)
rownames(tbl.j.rare) = tbl.j.rare$taxon
outFile = paste(c(opts[['-o']], 'rare'), collapse='_')
saveRDS(tbl.j.rare, outFile)
message('File written: ', outFile)
### dominant
tbl.j.dom = tbl.j %>% filter(mean_rel_abund_perc > abund.cut)
rownames(tbl.j.dom) = tbl.j.dom$taxon
outFile = paste(c(opts[['-o']], 'dom'), collapse='_')
saveRDS(tbl.j.dom, outFile)
message('File written: ', outFile)



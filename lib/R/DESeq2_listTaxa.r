#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: DESeq2_listTaxa.r [options] <DESeq2>

options:
  <DESeq2>         DESeq2 results object file ("-" if from STDIN).
  --padj=<q>       Max adjusted P-value (value < padj_cutoff).
  --padjBH=<qq>    Max adjusted P-value (BH method).
  --log2=<l>       Min log2FoldChange (log2fc).
  --inv            Get non-incorporators. 
  -h               Help
description:
  List taxa in a DESeq2 object.
  Use `padj`, `padjBD`, and/or `log2` to list
  just the incorporators. 
' -> doc

opts = docopt(doc)
padj.cut = as.numeric(opts[['--padj']])

# packages
pkgs <- c('dplyr', 'tidyr')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

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


## edit
### DESeq2 object
deseq.res = as.data.frame(deseq.res)
if(! 'taxon' %in% colnames(deseq.res)){
  deseq.res$taxon = rownames(deseq.res)
}

## writing taxa if no options provided
if(is.null(opts[['--padj']]) &&
   is.null(opts[['--padjBH']]) &&
   is.null(opts[['--log2']])){
  taxa = as.data.frame(deseq.res$taxon)
} else {
  ## calling incorporators
  ### log2fold cutoff & padj cutoff
  log2.cut = -1/0      # -negInfinity (ie., no cutoff)
  padj.cut = 1/0       # Infinity (ie., no cutoff)
  
  #### log2.cutoff
  if(! is.null(opts[['--log2neg']])){
    log2neg = as.numeric(opts[['--log2neg']])
    x = deseq.res %>%
      filter(log2FoldChange < 0) %>%
        mutate(log2FoldChange = abs(log2FoldChange)) %>%
          summarize(q5 = quantile(log2FoldChange, log2neg))
    log2.cut = abs(x[1,])
  }
  if(! is.null(opts[['--log2']])){
    log2.cut = as.numeric(opts[['--log2']])
  }
  message('Log2Fold cutoff: ', log2.cut)
  padj.cut = opts[['--padj']]
  #### padj.cut
  if(! is.null(opts[['--padj']])){
    padj.cut = as.numeric(opts[['--padj']])
    message('padj cutoff: ', padj.cut)
    deseq.res = deseq.res %>% 
      mutate(incorp = (padj < padj.cut) & (log2FoldChange >= log2.cut))
  } else if(! is.null(opts[['--padjBH']])){
    padj.cut = as.numeric(opts[['--padjBH']])
    message('padj cutoff: ', padj.cut)
    deseq.res = deseq.res %>% 
      mutate(incorp = (padj.BH < padj.cut) & (log2FoldChange >= log2.cut))
  } else {
    stop('Provide either --padj or --padjBH')
  }
  # filtering to just incorps
  deseq.res = deseq.res %>% filter(incorp != opts[['--inv']])
  taxa = as.data.frame(deseq.res$taxon)
}

con = pipe("cat", "wb")
write.table(taxa, con, row.names=FALSE, col.names=FALSE, quote=FALSE)

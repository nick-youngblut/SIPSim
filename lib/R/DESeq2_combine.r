#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: DESeq2_combine.r [options] <DESeq2>...

options:
  <DESeq2>         DESeq2 file(s)
  -o=<o>           Basename for output file.
                   [Default: DESeq2_comb]
  -h               Help
description:
  Combining >1 DESeq2 objects. P-values will be globally
  adjusted. 
' -> doc

opts = docopt(doc)


# packages
pkgs <- c('dplyr', 'tidyr')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# main
## import DESeq objects as a list
deseqs = sapply(opts[['<DESeq2>']], readRDS)
deseqs = lapply(deseqs, function(x) suppressPackageStartupMessages(as.data.frame(x)))

## combing objects
for (n in names(deseqs)){
  deseqs[[n]]$file = n
}
deseqs = do.call(rbind, deseqs)
deseqs$taxon = gsub('.+_DESeq2\\.', '', rownames(deseqs))
rownames(deseqs) = 1:nrow(deseqs)

## global adjustment of p-values
deseqs = deseqs %>%
  mutate(padj = p.adjust(p, 'BH')) %>%
    group_by(taxon) %>%
      mutate(log2FoldChange_Max = max(log2FoldChange)) %>%
        ungroup() %>%
          filter(log2FoldChange == log2FoldChange_Max) %>%
            select(-log2FoldChange_Max)

#deseqs %>% head %>% as.data.frame %>% print; stop()

## writing output
con = pipe("cat", "wb")
saveRDS(deseqs %>% as.data.frame, con)


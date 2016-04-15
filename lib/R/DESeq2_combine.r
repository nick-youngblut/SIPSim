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
  adjusted. For each taxon, the max log2 fold change (l2fc)
  of any DESeq2 object will be used. A "file" column
  will indicate which file was used for selecting l2fc.
' -> doc

opts = docopt(doc)

#print(padj_cut ); stop()

# packages
pkgs <- c('dplyr', 'tidyr')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# main
## import DESeq objects as a list
deseqs = sapply(opts[['<DESeq2>']], function(x) read.delim(x, sep='\t'))
deseqs = lapply(deseqs, function(x) suppressPackageStartupMessages(as.data.frame(x)))

## combing objects
for (n in names(deseqs)){
  deseqs[[n]]$file = n
}
deseqs = do.call(rbind, deseqs)
deseqs$taxon = gsub('.+_DESeq2\\.', '', rownames(deseqs))
rownames(deseqs) = 1:nrow(deseqs)


## global adjustment of p-values; using max log2FoldChange value
deseqs = deseqs %>%
  mutate(padj = p.adjust(p, 'BH')) %>%
    group_by(taxon) %>%
      mutate(log2FoldChange_Max = max(log2FoldChange)) %>%
        ungroup() %>%
          filter(log2FoldChange == log2FoldChange_Max) %>%
            select(-log2FoldChange_Max)

## writing output
con = pipe("cat", "wb")
saveRDS(deseqs %>% as.data.frame, con)


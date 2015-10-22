#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: comm_add_target.r [options] <comm> <target>

options:
  <comm>      Community file name ("-" if from stdin).
  <target>    Target taxa file name.
  -h          Help' -> doc

opts = docopt(doc)


# packages
pkgs <- c('dplyr', 'tidyr', 'phyloseq')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}



# main
## loading comm file
if(opts[['comm']] == '-'){
  df.comm = read.delim(file('stdin'), sep='\t')
} else {
  df.comm = read.delim(opts[['comm']], sep='\t')
}

## loading target file
df.target = read.delim(opts[['target']], sep='\t')


## adding targets
df.comm = left_join(df.comm, df.target, c('taxon_name' = 'OTU')) %>%
    as.data.frame %>%
    mutate(ssu_ID = ssu_ID %>% as.character,
           genome_fileID = genome_fileID %>% as.character,
           genomeID = genomeID %>% as.character,
           genome_seqID = genome_seqID %>% as.character,
           OTU_taxonomy = OTU_taxonomy %>% as.character)


# writing table
write.table(df.comm, stdout(), quote=FALSE, row.names=FALSE, sep='\t')

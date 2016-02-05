#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: comm_set_abund.r [options] <comm>

Options:
  <comm>      Community file name ("-" if from stdin).
  -s          Get abundances from subsampling existing taxa abundances.
  -h          Help

Description:
  Setting percent relative abundances to sum up to 100.
  In more technical terms: performing a closure operation.

  Ranks are also adjusted if needed (eg., if taxa were removed
  from the community table).
' -> doc

opts = docopt(doc)


# packages
pkgs <- c('dplyr', 'tidyr')
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


## closure operation
## rank adjusted (just to be sure)
df.comm = df.comm %>%
  group_by(library) %>%
    mutate(rel_abund_perc = rel_abund_perc / sum(rel_abund_perc) * 100,
           rank = row_number(-rel_abund_perc))



## writing table
write.table(df.comm, stdout(), quote=FALSE, row.names=FALSE, sep='\t')

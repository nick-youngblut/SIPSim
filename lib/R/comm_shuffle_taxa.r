#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: comm_shuffle_taxa.r [options] <comm>

Options:
  <comm>      Community file name ("-" if from stdin).
  -h          Help

Description:
  Shuffle taxa within each library of the community file.
  This can be used to assess how G+C content variation
  among the representative genomes affects the final HR-SIP
  simulated data.

  OUTPUT: community file written to STDOUT.

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


# shuffling taxa
df.comm = df.comm %>%
  group_by(library) %>%
    mutate(taxon_name = taxon_name %>% sample %>% sample)



## writing table
write.table(df.comm, stdout(), quote=FALSE, row.names=FALSE, sep='\t')

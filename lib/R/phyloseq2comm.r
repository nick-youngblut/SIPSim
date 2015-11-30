#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: phyloseq2comm.r [options] <phyloseq>

Options:
  <phyloseq>   Phyloseq object file.
  <outfile>    Output file name.
  -d=<d>       Select a day.
  -s=<s>       Select a substrate.
  -h           Help

Description:
  Converting a phyloseq object to a community
  table file
' -> doc

opts = docopt(doc)


# packages
pkgs <- c('dplyr', 'tidyr', 'phyloseq')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# functions
physeq2comm = function(physeq){
  df.otu = physeq %>% otu_table %>% as.matrix %>% as.data.frame
  df.otu$taxon_name = rownames(df.otu)
  df.otu = gather(df.otu, 'sample', 'abundance', 1:(ncol(df.otu)-1))
  df.otu = df.otu %>%
      group_by(sample) %>%
      mutate(rel_abund_perc = abundance / sum(abundance) * 100,
             rank = row_number(-abundance)) %>%
             ungroup() %>%
             mutate(library = gsub('^', '_', sample) %>% as.factor %>% as.numeric) %>%
             select(library, taxon_name, rel_abund_perc, rank) %>%
             arrange(library, rank)
  return(df.otu)
}


# main
## loading phyloseq object
physeq = readRDS(opts[['phyloseq']])
## pruning
day.bool = ! is.null(opts[['-d']])
sub.bool = ! is.null(opts[['-s']])
if(day.bool | sub.bool){
  physeq.m = physeq %>% sample_data
 if(day.bool & sub.bool){
   physeq = prune_samples(physeq.m$Substrate == opts[['-s']] &
                          physeq.m$Day == opts[['-d']], physeq) %>%
                filter_taxa(function(x) sum(x) > 0, TRUE)
 } else if(day.bool){
   physeq = prune_samples(physeq.m$Day == opts[['-d']], physeq) %>%
     filter_taxa(function(x) sum(x) > 0, TRUE)
 } else if(sub.bool){
   physeq = prune_samples(physeq.m$Substrate == opts[['-s']], physeq) %>%
     filter_taxa(function(x) sum(x) > 0, TRUE)
 }
}
## conversion
df = physeq2comm(physeq)
## writing table
write.table(df, stdout(), quote=FALSE, row.names=FALSE, sep='\t')


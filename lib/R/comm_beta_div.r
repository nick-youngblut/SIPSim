#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: comm_beta_div.r [options] <comm>

Options:
  <comm>    Community file name ("-" if from stdin).
  -b=<b>    Beta-diversity metrics. See vegan::vegdist for options.
            [Default: bray]
  -d        Write out diagonal distance matrix
            instead of column-formatted.
  -h        Help

Description:
  Calculate the beta-diversity among libraries (gradients).

  Output written to STDOUT.
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

# converting to a wide table
df.comm.w = df.comm %>%
  dplyr::select(library, taxon_name, rel_abund_perc) %>%
    spread(taxon_name, rel_abund_perc)
rownames(df.comm.w) = df.comm.w$library
df.comm.w$library = NULL

# calculating beta diversity
df.comm.beta = vegan::vegdist(df.comm.w, method=opts[['-b']]) 

## writing diagonal table
if(opts[['-d']]==TRUE){
  write.table(df.comm.beta %>% as.matrix %>% as.data.frame,
              stdout(), quote=FALSE, row.names=FALSE, sep='\t')
} else {
  # converting to column-format
  nam = df.comm.beta %>% as.matrix %>% as.data.frame %>% rownames
  df.col = combn(nam, 2) %>% t %>% as.data.frame
  df.col$beta_div = df.comm.beta %>% as.vector
  colnames(df.col) = c('library_x', 'library_y', opts[['-b']])
  
  # writing out table
  write.table(df.col, stdout(), quote=FALSE, row.names=FALSE, sep='\t')  
}

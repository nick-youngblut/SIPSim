#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: qSIP_confuseMtx.r [options] <BD_shift> <qSIP_atomExcess>

options:
  <BD_shift>         BD_shift file.
  <qSIP_atomExcess>  Output table from qSIP_atomExcess subcommand.
                     ("-" if from STDIN)
  -o=<o>             Basename for output files.
                     [Default: qSIP-cMtx]
  --BD=<b>           BD shift cutoff for identifying isotope incorporators.
                     [Default: 0.001]
  --libs=<l>         Libraries that are 13C treatments. (comma-separated list).
                     [Default: 2]
  -h                 Help
description: 
  Use caret to make a confusion matrix comparing
  KNOWN isotope incorporation (based on BD distribution
  shift between pre-incorp and post-incorp BD KDEs) and
  PREDICTED isotope incorporation (based on qSIP method).

  NOTE: all taxa with NA (or blank) in the `atom_CI_low` column,
  will be identified as non-incorporators.
' -> doc

opts = docopt(doc)
BD_shift.cut = as.numeric(opts[['--BD']])
libs = strsplit(opts[['--libs']], split=',')
libs = unlist(libs)

# packages
pkgs <- c('dplyr', 'tidyr', 'caret')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# main
## import
if(opts[['<qSIP_atomExcess>']] == '-'){
  con = pipe("cat", "rb")
  df_qSIP = read.delim(con, sep='\t')
} else {
  df_qSIP = read.delim(opts[['<qSIP_atomExcess>']], sep='\t')
}
df_shift = read.delim(opts[['BD_shift']], sep='\t')


## calling qSIP incorporators
### Note: if atom_CI_low is NA, then calling non-incorporator
### qSIP
df_qSIP = df_qSIP %>%
  mutate(incorp = ifelse(atom_CI_low > 0, TRUE, FALSE),
         incorp = ifelse(is.na(incorp), FALSE, incorp))


### BD-shift table (reference)
if (ncol(df_shift) == 8){
  df_shift = df_shift %>%
    mutate(incorp = median > BD_shift.cut)
} else {
  df_shift = df_shift %>%
    filter(lib2 == '2') %>%
      dplyr::rename('library' = lib2) %>%
        mutate(incorp = BD_shift > BD_shift.cut)
}
df_shift = df_shift %>%
  filter(library %in% libs)

# making factors of incorporation status
order_incorp = function(x){
  x$incorp = factor(x$incorp, levels=c(TRUE, FALSE))
  return(x)
}
df_shift = order_incorp(df_shift)
df_qSIP = order_incorp(df_qSIP)

# joining tables
df_shift$taxon = as.character(df_shift$taxon)
df_qSIP$taxon = as.character(df_qSIP$taxon)
df.j = inner_join(df_shift, df_qSIP, c('taxon' = 'taxon'))


# making confusion matrix
## incorp.x = df_shift  (KNOWN)
## incorp.y = qSIP      (PREDICTED)
df.j = df.j %>%
  mutate(incorp.known = incorp.x,
         incorp.pred = incorp.y) %>%
           select(-incorp.x, -incorp.y)

cfs.mtx = function(x){
  mtx = confusionMatrix(x$incorp.pred, x$incorp.known)
  return(mtx)
}

df.j.cmtx = df.j %>%
  group_by(library) %>%
    nest() %>%
      mutate(c.mtx = lapply(data, cfs.mtx)) 


get_tbl = function(obj, key){
  obj[[key]] %>% as.data.frame
}

conv = function(df){
  df = df %>% as.data.frame
  x = rownames(df) %>% as.data.frame
  df = cbind(x, df)
  colnames(df)[1] = 'variables'
  colnames(df)[2] = 'values'
  return(df)
}

df.j.table = df.j.cmtx %>%
  unnest(purrr::map(c.mtx, function(x) x[['table']] %>% as.data.frame)) %>%
    as.data.frame

df.j.overall = df.j.cmtx %>%
  unnest(purrr::map(c.mtx, function(x) x[['overall']] %>% conv)) %>%
    as.data.frame

df.j.byClass = df.j.cmtx %>%
  unnest(purrr::map(c.mtx, function(x) x[['byClass']] %>% conv)) %>%
    as.data.frame


## writing output
### raw data
saveRDS(df.j.cmtx, opts[['-o']])
message('File written: ', opts[['-o']])
outFile = paste(c(opts[['-o']], 'data.txt'), collapse='_')
write.table(df.j, outFile, sep='\t', quote=F, row.names=F)
message('File written: ', outFile)

## confusion matrix data
write_tbl = function(df, outPrefix, outFile){
  outFile = paste(c(outPrefix, outFile), collapse='_')
  write.table(df, file=outFile, sep='\t', quote=FALSE, row.names=FALSE)
  message('File written: ', outFile)
}
write_tbl(df.j.table, opts[['-o']], 'table.txt')
write_tbl(df.j.overall, opts[['-o']], 'overall.txt')
write_tbl(df.j.byClass, opts[['-o']], 'byClass.txt')


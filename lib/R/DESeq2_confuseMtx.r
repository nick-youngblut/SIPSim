#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: DESeq2_confuseMtx.r [options] <BD_shift> <DESeq2>

options:
  <BD_shift>       BD_shift file.
  <DESeq2>         DESeq2 results object file ("-" if from STDIN).
  -o=<o>           Basename for output files.
                   [Default: DESeq2-cMtx]
  --padj=<q>       Max adjusted P-value (value < padj_cutoff).
                   [Default: 0.1]
  --log2=<l>       Min log2FoldChange (log2fc). 
  --log2neg=<ln>   Min log2FoldChange based on negative log2fc values
                   (see below).
  --BD=<b>         BD shift cutoff for identifying isotope incorporators.
                   [Default: 0.001]
  --libs=<l>       Libraries that are 13C treatments. (comma-separated list).
                   [Default: 2]
  -h               Help
description:
  Use caret to make a confusion matrix comparing
  KNOWN isotope incorporation (based on BD distribution
  shift between pre-incorp and post-incorp BD KDEs) and
  PREDICTED isotope incorporation (based on DESeq2 results).

  log2neg
    This assumes that negative log2FoldChange (log2fc) values are caused
    by noise, and thus provide an estimate for the deviation
    in log2FoldChange from 0 that is the result of noise.
    The cutoff is determined as the quartile of negative log2fc values
    (converted to absolute values).
    For example: if log2neg=0.95, the log2fc cutoff is set as
    the distance from 0 that 95% of all negative "noise" log2fc values.    
' -> doc

opts = docopt(doc)
padj.cut = as.numeric(opts[['--padj']])
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
if(opts[['<DESeq2>']] == '-'){
  con = pipe("cat", "rb")
  #deseq.res = readRDS(con)
  deseq.res = read.delim(con, sep='\t')
} else {
  #deseq.res = readRDS(opts[['<DESeq2>']])
  deseq.res = read.delim(opts[['<DESeq2>']], sep='\t')
}
x = suppressPackageStartupMessages(as.data.frame(deseq.res))
BD.shift = read.delim(opts[['BD_shift']], sep='\t')

## edit
### DESeq2 object
deseq.res = as.data.frame(deseq.res)
if(! 'taxon' %in% colnames(deseq.res)){
  deseq.res$taxon = rownames(deseq.res)
}

### log2fold cutoff & padj cutoff
log2.cut = -1/0      # -negInfinity (ie., no cutoff)

#### log2.cutoff
message('Log2Fold cutoff: ', log2.cut)
message('padj cutoff: ', padj.cut)
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

#### padj.cut
deseq.res = deseq.res %>% 
  mutate(incorp = (padj < padj.cut) & (log2FoldChange >= log2.cut))



### BD-shift table (reference)
if (ncol(BD.shift) == 8){
  BD.shift = BD.shift %>%
    mutate(incorp = median > BD_shift.cut)
} else {
  BD.shift = BD.shift %>%
    filter(lib2 == '2') %>%
      dplyr::rename('library' = lib2) %>%
        mutate(incorp = BD_shift > BD_shift.cut)
}
BD.shift = BD.shift %>%
  filter(library %in% libs)

# making factors of incorporation status
order_incorp = function(x){
  x$incorp = factor(x$incorp, levels=c(TRUE, FALSE))
  return(x)
}
BD.shift = order_incorp(BD.shift)
deseq.res = order_incorp(deseq.res)

# joining tables
BD.shift$taxon = as.character(BD.shift$taxon)
deseq.res$taxon = as.character(deseq.res$taxon)
df.j = inner_join(BD.shift, deseq.res, c('taxon' = 'taxon'))


# making confusion matrix
## incorp.x = BD.shift  (KNOWN)
## incorp.y = deseq.res  (PREDICTED)
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



#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: qSIP_confuseMtx.r [options] <BD_shift> <qSIP_atomX>

Options:
  <BD_shift>       BD_shift file.
  <qSIP_atomX>     qSIP table with atom fraction excess confidence
                   intervalsv("-" if from STDIN).
  -o=<o>           Basename for output files.
                   [Default: qSIP-cMtx]
  -h               Help
Description:
  Use caret to make a confusion matrix comparing
  KNOWN isotope incorporation (based on BD distribution
  shift between pre-incorp and post-incorp BD KDEs) and
  PREDICTED isotope incorporation (based on qSIP method;
  Hungate et al., 2015)

References
  Hungate BA, Mau RL, Schwartz E, Caporaso JG, Dijkstra P, Gestel N van, et
  al. (2015). Quantitative Microbial Ecology Through Stable Isotope Probing.
  Appl Environ Microbiol AEM.02280-15
' -> doc

opts = docopt(doc)
#padj.cut = as.numeric(opts[['--padj']])

# packages
pkgs <- c('dplyr', 'caret')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# main
## import
if(opts[['<qSIP_atomX>']] == '-'){
  con = pipe("cat", "rb")
  deseq.res = readRDS(con)
} else {
  deseq.res = readRDS(opts[['<qSIP_atomX>']])
}
x = suppressPackageStartupMessages(as.data.frame(deseq.res))
BD.shift = read.delim(opts[['BD_shift']], sep='\t')


#--------------
# TODO: finish
#--------------

## edit
### DESeq2 object
deseq.res = as.data.frame(deseq.res)
deseq.res$taxon = rownames(deseq.res)

#deseq.res %>% head %>% print; stop()

### log2fold cutoff & padj cutoff
log2.cut = -1/0   # -negInfinity (ie., no cutoff)
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
} else
if(! is.null(opts[['--padjBH']])){
  padj.cut = as.numeric(opts[['--padjBH']])
  message('padj cutoff: ', padj.cut)
  deseq.res = deseq.res %>% 
    mutate(incorp = (padj.BH < padj.cut) & (log2FoldChange >= log2.cut))
} else {
  stop('Provide either --padj or --padjBH')
}

### BD-shift table (reference)
BD.shift = BD.shift %>%
  filter(lib2 == '2') %>%
      mutate(incorp = BD_shift > 0.05)

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
tbl.c = inner_join(BD.shift, deseq.res, c('taxon' = 'taxon'))

# making confusion matrix
## incorp.x = BD.shift  (KNOWN)
## incorp.y = deseq.res  (PREDICTED)
tbl.c = tbl.c %>%
  mutate(incorp.known = incorp.x,
         incorp.pred = incorp.y) %>%
           select(-incorp.x, -incorp.y)
#    rename(c('incorp.x' = 'incorp.known', 'incorp.y' = 'incorp.pred'))
c.mtx = confusionMatrix(tbl.c$incorp.pred, tbl.c$incorp.known)

## writing
outFile = paste(c(opts[['-o']], 'data.csv'), collapse='_')
write.csv(tbl.c, outFile, quote=F, row.names=F)
message('File written: ', outFile)
saveRDS(c.mtx, opts[['-o']])
message('File written: ', opts[['-o']])

write_tbl = function(obj, key, outFile){
  x = obj[[key]] %>% as.data.frame()
  colnames(x)[1] = key
  write.csv(x, file=outFile, quote=F)
  message('File written: ', outFile)
}
outFile = paste(c(opts[['-o']], 'table.csv'), collapse='_')
write_tbl(c.mtx, 'table', outFile)
outFile = paste(c(opts[['-o']], 'overall.csv'), collapse='_')
write_tbl(c.mtx, 'overall', outFile)
outFile = paste(c(opts[['-o']], 'byClass.csv'), collapse='_')
write_tbl(c.mtx, 'byClass', outFile)


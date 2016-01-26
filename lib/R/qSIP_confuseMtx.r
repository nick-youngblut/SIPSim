#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: qSIP_confuseMtx.r [options] <BD_shift> <qSIP_atomX> <exp_design>

Options:
  <BD_shift>     BD_shift file.
  <qSIP_atomX>   qSIP table with atom fraction excess confidence
                 intervals ("-" if from STDIN).
  <exp_design>   Experimental design file (see qSIP_atomExcess).
  --CI=<c>       Confidence interval cutoff for identifying
                 isotope incorporators.
                 [Default: 0]
  --BD=<b>       BD shift cutoff for identifying isotope
                 incorporators.
                 [Default: 0.05]
  -o=<o>         Basename for output files.
                 [Default: qSIP-cMtx]
  -h             Help
Description:
  Use caret to make a confusion matrix comparing KNOWN isotope
  incorporation (based on BD distribution shift between pre-incorp
  and post-incorp BD KDEs) and PREDICTED isotope incorporation
  (based on qSIP method; Hungate et al., 2015).

References
  Hungate BA, Mau RL, Schwartz E, Caporaso JG, Dijkstra P, Gestel N van, et
  al. (2015). Quantitative Microbial Ecology Through Stable Isotope Probing.
  Appl Environ Microbiol AEM.02280-15
' -> doc

opts = docopt(doc)
CI.cut = as.numeric(opts[['--CI']])
BD_shift.cut = as.numeric(opts[['--BD']])


# packages
pkgs <- c('dplyr', 'caret')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# main

## import
### atom excess
if(opts[['<qSIP_atomX>']] == '-'){
  con = pipe("cat", "rb")
  df.atom = read.delim(con, sep='\t')
} else {
  df.atom = read.delim(opts[['<qSIP_atomX>']], sep='\t')
}
### BD shift
df.shift = read.delim(opts[['<BD_shift>']], sep='\t')
### exp. design
df.design = read.delim(opts[['<exp_design>']], sep='\t', header=FALSE)


##summarizing BD shift for treatments
df.design$V1 = sapply(df.design$V1, tolower)
treatments = df.design %>%
  filter(V2 == 'treatment')
treatments = treatments$V1

df.shift.s = df.shift %>%
  filter(lib2 %in% treatments) %>%
    group_by(taxon) %>%
      summarize(mean_BD_shift = mean(BD_shift, na.rm=TRUE))

## table join
df.j = inner_join(df.atom, df.shift.s, c('taxon' = 'taxon'))


## identifying incorporators
df.j = df.j %>%
  mutate(qSIP_incorper = ifelse(atom_CI_low > CI.cut, TRUE, FALSE),
         hrSIP_incorper = ifelse(mean_BD_shift > BD_shift.cut, TRUE, FALSE))


## confusion matrix
c.mtx = confusionMatrix(df.j$qSIP_incorper, df.j$hrSIP_incorper)


## writing
outFile = paste(c(opts[['-o']], 'data.csv'), collapse='_')
write.csv(df.j, outFile, quote=F, row.names=F)
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

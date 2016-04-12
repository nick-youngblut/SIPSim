#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: phyloseq_DESeq2.r [options] <phyloseq>

options:
  <phyloseq>     Phyloseq object file.
  --log2=<l>     Log2 fold change cutoff.
                 [Default: 0.25]
  --hypo=<h>     Alternative hypothesis tested by DESeq
                 ("greaterAbs","greater","less")
                 [Default: NULL]
  --cont=<c>     Control libraries. (comma-separated list)
                 [Default: 1]
  --treat=<t>    Treatment libraries. (comma-separated list)
                 [Default: 2]
  --label=<x>    Add a label to the DESeq2 object, which is
                 useful if combining multiple DESeq2 objects.
                 [Default: NULL]
  --occur=<oc>   Minimum fraction of samples that a taxon must
                 be found (ie., sparcity threshold).
                 If a comma-separated list is provided, each
                 cutoff is tested and the cutoff with the most
                 rejected hypotheses is used.
                 [Default: 0]
  --padj=<p>     Adjusted P-value cutoff for determining rejected
                 hypotheses (see `occur` option).
                 [Default: 0.1]
  --nrej=<nr>    Output file name for table with number of rejected
                 hypotheses per `occur` cutoff.
                 NOTE: DESeq2 output written to STDOUT.
                 [Default: DESeq2_n-rejectHypos.txt]
  -h             Help

description:
  Run DESeq2 on a phyloseq file produced in the SIPSim pipeline.
  DESeq2 output written to STDOUT.
' -> doc

opts = docopt(doc)
cont = unlist(strsplit(opts[['--cont']], split=','))
treat = unlist(strsplit(opts[['--treat']], split=','))
occurs = unlist(strsplit(opts[['--occur']], split=','))
occurs = as.vector(sapply(occurs, as.numeric))
padj.cut = as.numeric(opts[['--padj']])
log2.cut = as.numeric(opts[['--log2']])
nrej_out_file = opts[['--nrej']]

# packages
pkgs <- c('phyloseq', 'DESeq2')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}


#-- functions --#
prune_by_occur = function(physeq, occur){
  # taxa must be in `occur` number of gradient fractions
  if(is.null(occur)){
    occur = 0
  }
  msg = paste0(' pre-filter: number of taxa:', ntaxa(physeq))
  write(msg, stderr())
  physeq.prn = filter_taxa(physeq, function(x) sum(x > 0) > (occur * length(x)), TRUE)
  msg = paste0(' post-filter: number of taxa:', ntaxa(physeq.prn))
  write(msg, stderr())
  return(physeq.prn)
}

gm_mean = function(x, na.rm=TRUE){
  # calculate the geometric mean
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }


make_condition = function(dds){
  # making 'condition' factor for DESeq2
  condition = as.character(dds$condition)
  for (x in cont){
    condition = ifelse(condition == x, 'control', condition)
  }
  for (x in treat){
    condition = ifelse(condition == x, 'treatment', condition)
  }

  dds$condition = factor(condition, levels=c('control', 'treatment'))
  
  if(any(is.na(dds$condition))){
    stop('\nSome control/treatments are NA. Check --cont & --treat.\n')
  }
  return(condition)
}


deseq_run = function(occur, physeq, opts){
  # deseq analysis on samples
  ## status
  msg = paste0('\n### DESeq2 with occurance cutoff: ', occur)
  write(msg, stderr())  
  
  ## prune physeq by occurance
  physeq.prn = prune_by_occur(physeq, occur)

  ## making deseq2 object
  dds = phyloseq_to_deseq2(physeq.prn, ~condition)
  physeq.prn = NA

  ## making/setting 'condition' factor
  condition = make_condition(dds)
  
  ## calculate geometric means prior to estimate size factors
  ### This method is not sensitive to zeros
  geoMeans = apply(counts(dds), 1, gm_mean)
  dds = estimateSizeFactors(dds, geoMeans=geoMeans)

  ## DESeq2 call
  dds = DESeq(dds, fitType='local', test='Wald')
  if(opts[['--hypo']] != 'NULL'){
    res = results(dds,
      lfcThreshold=log2.cut,
      altHypothesis=opts[['--hypo']],
      independentFiltering=FALSE)  
  } else {
    res = results(dds, independentFiltering=FALSE)
  }

  ## p-adjust as in Pepe-Ranney et al. (2015)
  beta = res$log2FoldChange
  betaSE = res$lfcSE
  res$p = pnorm(beta, log2.cut, betaSE, lower.tail = FALSE)
  res$padj.BH = p.adjust(res$p, "BH")
  res$label = opts[['--label']]

  ## return
  return(res)
}

n_reject = function(dds_res, padj.cut){
  # finding the number of hypotheses rejected
  length(dds_res$padj[dds_res$padj < padj.cut])
}


#-- main --#
# import
if(opts[['<phyloseq>']] == '-'){
  con = pipe("cat", "rb")
  physeq = readRDS(con)
} else {
  physeq = readRDS(opts[['<phyloseq>']])
}

# library must be a character
physeq.sd = sample_data(physeq)
physeq.sd$library = as.factor(as.character(physeq.sd$library))
physeq.sd$condition = physeq.sd$library
sample_data(physeq) = physeq.sd

# DESeq2 runs for each occurance (sparsity) cutoff
## returns a list of DESeq2 results objects
res_all = sapply(occurs, deseq_run, physeq=physeq, opts=opts)

# Finding number of hypotheses rejected
n_rej = lapply(res_all, n_reject, padj.cut=padj.cut)
n_rej = unlist(n_rej)
df_n_rej = data.frame(occur = occurs, n_hypo_rejected = n_rej)
write.table(df_n_rej, nrej_out_file, sep='\t', quote=FALSE, row.names=FALSE)
msg = paste0('\n# N-rejected hypothesis table written: ', nrej_out_file)
write(msg, stderr())


# selecting DESeq2 results with most number of rejected Hs
max_n_rej = max(n_rej)
max_n_rej_i = which(n_rej == max_n_rej)
if(length(max_n_rej_i) > 1){
  max_n_rej_i = max_n_rej_i[1]
}
res = res_all[[max_n_rej_i]]
## status
occur_best = occurs[max_n_rej_i]
msg = paste0('# `occur` cutoff with max number of rejected hypotheses: ',
  occur_best)
write(msg, stderr())
msg = paste0('# The max number of rejected hypotheses: ',
  max_n_rej, '\n')
write(msg, stderr())


## writing
con = pipe("cat", "wb")
saveRDS(res, con)

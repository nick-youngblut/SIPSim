#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: comm_add_richness.r [options] <comm> <richness>

Options:
  <comm>      Community file name ("-" if from stdin).
  <richness>  A table of richness for each community (library).
  -s          Get abundances from subsampling existing taxa abundances.
  -h          Help

Description:
  Adding random OTUs (taxa) to community table
  to reach the richness as specified in the <richness> table.

  Abundances of new "random" taxa: 
    Abundances of new random taxa are obtained by either:
      i) curve-fitting on the abundances in the provided comm
         file (rel_abund_perc)
      ii) if "-s": selecting abundances from subsampling
          (with replacement) of the abundances in the comm file

  The random taxa will be appended to the end of the provided
  community table and written to STDOUT.

  richness table:
    2-column tab-delimited table: library_id(TAB)richness
    No header!
    "richness" indicates total community richness needed.
' -> doc

opts = docopt(doc)


# packages
pkgs <- c('dplyr', 'tidyr', 'fitdistrplus')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# functions
fit_dist = function(x){
  # fitting a distribution to the taxon relative abundances
  x = as.numeric(as.vector(x))
  # fit
  f.norm = fitdist(x, 'norm')
  f.exp = fitdist(x, 'exp')
  f.logn = fitdist(x, 'lnorm')
  f.gamma = fitdist(x, 'gamma')
  
  f.list = list(f.norm, f.exp, f.logn, f.gamma)
  names = c('normal', 'exponential', 'lognormal', 'gamma')
  # best fit (via AIC)
  stat.res = gofstat(f.list, fitnames=names)
  best.AIC = which(stat.res$aic == min(stat.res$aic))
  best.AIC.func = f.list[[best.AIC[1]]]
  
  if (length(best.AIC.func$estimate) > 1){
    best.AIC.func$estimate = best.AIC.func$estimate[1]
  }
  msg = paste0('Best curve fit of abundances: ', names[best.AIC[1]])
  write(msg, stderr())
  msg = paste0('  Estimate: ', best.AIC.func$estimate)
  write(msg, stderr())
  return(best.AIC.func)
}


get_rel_abunds = function(abund.dist.func, richness){
  # sampling from the fitted abundance distribution
  func.l = list(
    'norm' = rnorm,
    'exp' = rexp,
    'lnorm' = rlnorm,
    'gamma' = rgamma)

  func.name = abund.dist.func$distname
  func = func.l[[func.name]]
  abunds = vector()
  if (func.name == 'exponential'){  # 1 parameter function
    abunds = func(richness, abund.dist.func$estimate)
  } else {
    abunds = func(richness, abund.dist.func$estimate, abund.dist.func$sd)
  }
  return(abunds)
}

make_rand_comm = function(comm.cols, abund.dist.func, lib, richness){
  # getting random taxa abundances from abundance distribution
  ## matrix of random taxa
  m = matrix(ncol=length(comm.cols),
    nrow=richness,
    dimnames=list(c(), comm.cols))
  ## lib
  m[,'library'] = lib
  ## OTU IDs
  m[,'taxon_name'] = gsub('^', 'rand.', 1:richness)
  ## relative abundances
  m[,'rel_abund_perc'] = get_rel_abunds(abund.dist.func, richness)
  m = as.data.frame(m)
  m$rel_abund_perc = as.numeric(as.vector(m$rel_abund_perc))
  return(m)
}

make_rand_comm_subsample = function(comm.cols, abunds,  lib, richness){
  # getting random taxa abundances from subsampling of other taxa
  m = matrix(ncol=length(comm.cols),
    nrow=richness,
    dimnames=list(c(), comm.cols))
    
  ## lib
  m[,'library'] = lib
  ## OTU IDs
  m[,'taxon_name'] = gsub('^', 'rand.', 1:richness)
  ## relative abundances
  m[,'rel_abund_perc'] = sample(abunds, richness, replace=TRUE)
  m = as.data.frame(m)
  m$rel_abund_perc = as.numeric(as.vector(m$rel_abund_perc))
  return(m)
}


# main
## loading comm file
if(opts[['comm']] == '-'){
  df.comm = read.delim(file('stdin'), sep='\t')
} else {
  df.comm = read.delim(opts[['comm']], sep='\t')
}

## loading richness table
df.rich = read.delim(opts[['richness']], sep='\t', header=FALSE)


## for each library, adding to richness
df.comm.cols = colnames(df.comm)
for (i in 1:nrow(df.rich)){
  lib = as.character(df.rich[i,1])
  need_richness = df.rich[i,2]

  # current community richness
  curr_richness = length(unique(df.comm[df.comm[,'library'] == lib,'taxon_name']))
  diff_richness = need_richness - curr_richness
  
  # status
  msg = paste0('Library: ', lib,
    '\n  Current richness: ', curr_richness,
    '\n  Richness needed: ', need_richness,
    '\n  Additional taxa needed: ', diff_richness)
  write(msg, stderr())

  if (diff_richness <= 0){
    msg = paste0('    Additional taxa needed <= 0. Skipping')
    write(msg, stderr())
    next
  }
  
  # curve fitting relative abundance
  if (opts[['-s']] == TRUE){
    abunds = df.comm[df.comm[,'library'] == lib,'rel_abund_perc']
    df.rand = make_rand_comm_subsample(df.comm.cols, abunds, lib, diff_richness)
    df.comm = rbind(df.comm, df.rand)
  } else {
    df.comm.p = filter(df.comm, library == lib)
    abund.dist.func = fit_dist(df.comm.p$rel_abund_perc)

    # richness accounting for existing taxa
    df.rand = make_rand_comm(df.comm.cols, abund.dist.func, lib, diff_richness)
    df.comm = rbind(df.comm, df.rand)      
  }
  ## status
  msg = paste0('  CHECK: final richness: ', length(df.comm$taxon_name))
  write(msg, stderr())
}

## adjusting rank
df.comm = df.comm %>%
  group_by(library) %>%
    mutate(rank = row_number(-rel_abund_perc))



## writing table
write.table(df.comm, stdout(), quote=FALSE, row.names=FALSE, sep='\t')

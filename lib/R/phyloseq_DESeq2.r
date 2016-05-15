#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: phyloseq_DESeq2.r [options] <phyloseq>

options:
  <phyloseq>          Phyloseq object file.
  -w=<w>              "heavy" BD window(s). (omma separated list)
                      [Default: 1.71-1.75]
  --log2=<l>          Log2 fold change cutoff.
                      [Default: 0.25]
  --hypo=<h>          Alternative hypothesis tested by DESeq
                      ("greaterAbs","greater","less")
                      [Default: greater]
  --cont=<c>          Control libraries. (comma-separated list)
                      [Default: 1]
  --treat=<t>         Treatment libraries. (comma-separated list)
                      [Default: 2]
  --occur_all=<oa>    Minimum fraction of samples that a taxon must
                      be found (ie., sparcity threshold).
                      If a comma-separated list is provided, each
                      cutoff is tested and the cutoff with the most
                      rejected hypotheses is used. Cutoff applied
                      to ALL samples.
                      [Default: 0]
  --occur_heavy=<oh>  Minimum fraction of samples that a taxon must
                      be found (ie., sparcity threshold).
                      If a comma-separated list is provided, each
                      cutoff is tested and the cutoff with the most
                      rejected hypotheses is used. Cutoff applied
                      to JUST HEAVY samples.
                      [Default: 0]
  --padj=<p>          Adjusted P-value cutoff for determining rejected
                      hypotheses (see the `occur` options).
                      [Default: 0.1]
  --all=<a>           Output file name for table with all DESeq2 results
                      for all DESeq2 runs.
                      [Default: DESeq2_all_runs.txt]
  -h                  Help

description:
  Run DESeq2 on a phyloseq file produced in the SIPSim pipeline.

  ALGORITHM:
  * All pairwise runs of DESeq2 for each occur/BD_window combination.
    * log2 fold change for each taxon between control & treatment
  * For each BD_window, select the occur cutoff with the greatest
    number of rejected hypotheses.
    * If multiple have the same max rej-hypos, then the lowest is used.
  * Globally adjust p-values for all multiple hypotheses across
    all BD windows (windows with the "Best" occur cutoff).
    * eg., 3 BD windows means 3 tests per taxon.
    * Benjamini-Hochberg correction
  * For each taxon, select the BD window with the highest
    log2 fold change value.

  An occurance (sparsity) cutoff can be applied before and/or
  after pruning samples to just "heavy" gradient fraction samples.

  CONT,TREAT:
    To use all libraries for the control and treatment, provide them all.
    If not all libraries are provided, the phyloseq object will be
    subsetted to just the user-provided libraries.

  OUTPUT:
    A table of all DESeq2 results is written to a file (see `all` option).
    A DESeq2 results table of just the "best" occurance cutoff and BD
    window log2 fold change values is written to STDOUT.

' -> doc

opts = docopt(doc)
opts[['--cont']] = unlist(strsplit(opts[['--cont']], split=','))
opts[['--treat']] = unlist(strsplit(opts[['--treat']], split=','))
opts[['--padj']] = as.numeric(opts[['--padj']])
opts[['--log2']] = as.numeric(opts[['--log2']])

# occurance (sparsity)
opts[['--occur_all']] = unlist(strsplit(opts[['--occur_all']], split=','))
opts[['--occur_heavy']] = unlist(strsplit(opts[['--occur_heavy']], split=','))

# BD windows
opts[['-w']] = unlist(strsplit(opts[['-w']], split=','))


# packages
pkgs <- c('phyloseq', 'DESeq2', 'dplyr')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}


#-- functions --#
gm_mean = function(x, na.rm=TRUE){
  # calculate the geometric mean
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }


set_condition = function(dds, opts){
  cont = opts[['--cont']]
  treat = opts[['--treat']]
  
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
  return(dds)
}

prune_by_BD_window = function(BD_min, BD_max, physeq){
  # prune physeq object to just samples in BD window
  # status
  msg = paste0('# "heavy" BD window: ', BD_min, '-', BD_max)
  write(msg, stderr())  
  msg = paste0(' pre-filter: number of samples: ', nsamples(physeq))
  write(msg, stderr())

  # pruning to heavy window samples
  physeq.sd = sample_data(physeq)
  physeq.prn = prune_samples((physeq.sd$BD_mid >= BD_min) &
                               (physeq.sd$BD_mid <= BD_max), physeq)  

  # status
  msg = paste0(' post-filter: number of samples: ', nsamples(physeq.prn))
  write(msg, stderr())
  
  return(physeq.prn)
}


prune_by_occur = function(occur_cut, physeq){
  # pruning taxa by occurance
  msg = paste0('# occurance cutoff: ', occur_cut)
  write(msg, stderr())  

  msg = paste0(' pre-filter: number of taxa: ', ntaxa(physeq))
  write(msg, stderr())
  
  func = function(x) sum(x > 0) > (occur_cut * length(x))
  physeq.prn = filter_taxa(physeq, func, TRUE)

  msg = paste0(' post-filter: number of taxa: ', ntaxa(physeq.prn))
  write(msg, stderr())
  
  return(physeq.prn)
}


deseq_run = function(params, physeq, opts){
  # DESeq2 call on phyloseq object
  ## unpacking parameters
  occur_all = as.numeric(params[1])
  occur_heavy = as.numeric(params[2])
  BD_min = as.numeric(params[3])
  BD_max = as.numeric(params[4])
  
  # prune taxa by occurance cutoff (all samples)
  physeq.prn = prune_by_occur(occur_all, physeq)

  # prune sample to just BD window
  physeq.prn = prune_by_BD_window(BD_min, BD_max, physeq.prn)

  # prune taxa by occurance cutoff (heavy window samples)
  physeq.prn = prune_by_occur(occur_heavy, physeq.prn)
      
  # deseq analysis on samples
  ## status
  write('# DESeq2 run', stderr())  
  
  ## making deseq2 object
  dds = phyloseq_to_deseq2(physeq.prn, ~condition)
  physeq.prn = NA

  ## making/setting 'condition' factor
  dds = set_condition(dds, opts)
  if(dds$condition %>% unique %>% length == 1){
    msg = 'WARNING: no taxa present for either control or treatment. Skipping'
    write(msg, stderr())  
    return(NA)
  }
  
  ## calculate geometric means prior to estimate size factors
  ### This method is not sensitive to zeros
  geoMeans = apply(counts(dds), 1, gm_mean)
  dds = estimateSizeFactors(dds, geoMeans=geoMeans)

  ## DESeq2 call
  dds = DESeq(dds, fitType='local')
  res = results(dds, independentFiltering=FALSE)

  ## p-adjust as in Pepe-Ranney et al. (2015)
  beta = res$log2FoldChange
  betaSE = res$lfcSE
  res$p = pnorm(beta, opts[['--log2']], betaSE, lower.tail=FALSE)
  res$padj = p.adjust(res$p, "BH")
  ## setting all NA padj to 1 
  res$padj = ifelse(is.na(res$padj), 1, res$padj)
  
  ## adding pruning info to object
  res$occur_all = occur_all
  res$occur_heavy = occur_heavy
  res$heavy_BD_min = BD_min
  res$heavy_BD_max = BD_max
  res$taxon = rownames(res)
  rownames(res) = 1:nrow(res)

  # providing status on number of rejected hypothese
  n_rej_hypo = sum(res$padj < opts[['--padj']], na.rm=TRUE)
  msg = paste0('# Number of rejected hypotheses: ', n_rej_hypo)
  write(msg, stderr())  
  
  ## return
  write('', stderr())
  return(res)
}


write_maxRej_summary = function(deseq_res, opts){
  # writing a summary table of max rejected hypotheses for each BD window
  deseq_res_max = deseq_res %>%
    group_by(heavy_BD_min, heavy_BD_max, occur_all, occur_heavy) %>%    
      summarize(n_rej_hypo = sum(padj < opts[['--padj']], na.rm=TRUE)) %>%
        group_by(heavy_BD_min, heavy_BD_max) %>%
          summarize(max_n_rej_hypo = max(n_rej_hypo),
                    min_occur_all = min(occur_all),
                    max_occur_all = max(occur_all),
                    min_occur_heavy = min(occur_heavy),
                    max_occur_heavy = max(occur_heavy)) %>%
                      ungroup %>%
                        as.data.frame

  if(any(deseq_res_max$min_occur_all != deseq_res_max$min_occur_all)){
    stop('Occur_all cutoff logic error!')
  }
  if(any(deseq_res_max$min_occur_heavy != deseq_res_max$min_occur_heavy)){
    stop('Occur_heavy cutoff logic error!')
  }
  deseq_res_max = deseq_res_max %>% select(-max_occur_all, -max_occur_heavy)
  
  msg = paste0('\n#- Greatest number of rejected hypotheses for each BD window -#')
  write(msg, stderr())
  msg = paste0(c('heavy_BD_min', 'heavy_BD_max', 'max_n_rej_hypo', 'occur_all', 'occur_heavy'),
    collapse='\t')
  write(msg, stderr())
  df_to_stderr = function(x) write(paste0(x, collapse='\t'), stderr())
  blnk = apply(deseq_res_max, 1, df_to_stderr)
  write('', stderr())
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

# filtering to user-porvided control/treatment libraries
physeq.libs = physeq.sd$library %>% unique %>% as.vector
user.libs = c(opts[['--cont']], opts[['--treat']])
user.libs.only = setdiff(user.libs, physeq.libs)
if(length(user.libs.only) > 0){
  y = paste(user.libs.only, collapse=',')
  x = paste(c('Cannot find libraries:',y,'in the phyloseq object'), collapse=' ')
  stop(x)
}
physeq.libs.only = setdiff(physeq.libs, user.libs)
if(length(physeq.libs.only) > 0){
  y = paste(intersect(user.libs, physeq.libs), collapse=',')
  msg = paste(c('NOTE: filtering phyloseq to libraries: ', y, '\n'), collapse=' ')
  write(msg, stderr())
  physeq = subset_samples(physeq, library %in% user.libs)
}

# all pairwise of occur_all, BD_window, post-occur
var_params = expand.grid(opts[['--occur_all']],
  opts[['-w']],
  opts[['--occur_heavy']])
colnames(var_params) = c('occur_all', 'BD_window', 'occur_heavy')
var_params$BD_min = gsub('-.+', '', var_params$BD_window) %>% as.numeric
var_params$BD_max = gsub('.+-', '', var_params$BD_window) %>% as.numeric
var_params$BD_window = NULL

# DESeq2 call
deseq_res = apply(var_params, 1, deseq_run, physeq=physeq, opts=opts)
deseq_res = deseq_res[!is.na(deseq_res)]
deseq_res = do.call(rbind, deseq_res) %>%
  as.data.frame %>%
    filter(!is.na(taxon))

# write out full table of all DESeq2 data
write.table(deseq_res, opts[['--all']], sep='\t', quote=FALSE, row.names=FALSE)
msg = paste0('# A file of all DESeq2 run results was written to: ', opts[['--all']])
write(msg, stderr())


# Finding number of hypotheses rejected; selecting occurances with most rejects
## occur cutoffs with max number of rejected hypos
## if multiple cutoffs with the same rejects, then selecting the first
deseq_res = deseq_res %>%
  group_by(occur_all, occur_heavy, heavy_BD_min, heavy_BD_max) %>%
    mutate(n_rej_hypo = sum(padj < opts[['--padj']], na.rm=TRUE)) %>%
      group_by(heavy_BD_min, heavy_BD_max) %>%
        filter(n_rej_hypo == max(n_rej_hypo)) %>%
          ungroup() %>%
            distinct(heavy_BD_min, heavy_BD_max, taxon)

# writing table of max number of rejected hypos per BD window
write_maxRej_summary(deseq_res, opts)

# global P-value adjustment for all BD windows
msg = '# Global adjustment of p-values for all BW windows with "best" occurance cutoff'
write(msg, stderr())
deseq_res = deseq_res %>%
  mutate(padj_global = p.adjust(deseq_res$pvalue, method = "BH"),
         padj = ifelse(padj_global >= padj, padj_global, padj)) %>%
           select(-padj_global)

# writing table of max number of rejected hypos per BD window
## post global p-value adjustment
write_maxRej_summary(deseq_res, opts)

## for each taxon, selecting heavy window with the highest l2fc
deseq_res = deseq_res %>%
  group_by(taxon) %>%
    filter(log2FoldChange == max(log2FoldChange, na.rm=TRUE)) %>%
      ungroup() %>%
        distinct(taxon) %>%
          select(-n_rej_hypo) %>%
            as.data.frame

## writing table of results
con = pipe("cat", "wb")
write.table(deseq_res, con, sep='\t', quote=FALSE, row.names=FALSE)


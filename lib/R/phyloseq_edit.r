#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: phyloseq_edit.r [options] <phyloseq>

options:
  <phyloseq>         Phyloseq object file.
  --sampleSum=<ss>   Minimum sample size.
                     [Default: 100]
  --occur=<oc>       Minimum fraction of samples that a taxon must
                     be found (ie., sparcity threshold).
                     [Default: 0]
  --BD_min=<min>     Minimum sample buoyant density.
                     [Default: 1]
  --BD_max=<max>     Maximum sample buoyant density.
                     [Default: 2]
  --taxa=<t>         A file listing taxa to keep (and remove
                     the rest). One taxon per line.
  -h                 Help

description:
  Filter a phyloseq object to just samples/taxa of interest.

  Output: written to STDOUT.
' -> doc

opts = docopt(doc)


# packages
pkgs <- c('phyloseq')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# main
## import
if(opts[['<phyloseq>']] == '-'){
  con = pipe("cat", "rb")
  physeq = readRDS(con)
} else {
  physeq = readRDS(opts[['<phyloseq>']])
}

## set variables
sampleSum = as.numeric(opts[['--sampleSum']])
occur = as.numeric(opts[['--occur']])
BD_min = as.numeric(opts[['--BD_min']])
BD_max = as.numeric(opts[['--BD_max']])
taxa_file = opts[['--taxa']]


## Prune to select taxa (if needed
if(! is.null(taxa_file)){
  write('Pruning taxa...', stderr())
  taxa_to_keep = read.table(taxa_file, sep='\t')
  taxa_to_keep = as.character(taxa_to_keep[,1])
  preFilt_ntaxa = ntaxa(physeq)
  physeq = prune_taxa(taxa_to_keep, physeq)
  postFile_ntaxa = ntaxa(physeq)
  msg = paste0(' pre-filter: number of taxa:', preFilt_ntaxa)
  write(msg, stderr())
  msg = paste0(' post-filter: number of taxa:', postFile_ntaxa)
  write(msg, stderr())
}


## Pruning samples
write('Pruning samples...', stderr())
msg = paste0(' pre-filter: number of samples:', nsamples(physeq))
write(msg, stderr())
msg = paste0(' pre-filter: number of taxa:', ntaxa(physeq))
write(msg, stderr())

### sampleSum
if(! is.null(sampleSum)){
  physeq = prune_samples(sample_sums(physeq) >= sampleSum, physeq)
}
### occurance
if(! is.null(occur)){
  physeq = filter_taxa(physeq, function(x) sum(x > 0) > (occur * length(x)), TRUE)
}

### BD min-max
if(! is.null(BD_min) & ! is.null(BD_max)){
  physeq.sd = sample_data(physeq)
  physeq = prune_samples((physeq.sd$BD_mid >= BD_min) &
                           (physeq.sd$BD_mid <= BD_max), physeq)                            
  }

# status
msg = paste0(' post-filter: number of samples:', nsamples(physeq))
write(msg, stderr())
msg = paste0(' post-filter: number of taxa:', ntaxa(physeq))
write(msg, stderr())


## writing
con = pipe("cat", "wb")
saveRDS(physeq, con)

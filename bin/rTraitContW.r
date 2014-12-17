#!/usr/bin/Rscript
mod <- "5/19/14";
version <- "0.1";
author <- "Nick Youngblut";

### packages
pkgs <- c('optparse','ape')
for(i in 1:length(pkgs)){
  tmp <- pkgs[i]
  suppressPackageStartupMessages(library(pkgs[i], character.only=TRUE))
}

### I/O
# initialize options
option_list <- list(
  make_option(c("-t", "--tree"), type="character", help="Tree file"),
  make_option(c("-f", "--format"), type="character", default="newick", help="Tree file format (newick or nexus). [newick]"),
  make_option(c("-w", "--weight"), type="numeric", default=1, help="weight parameter (should be from 0-1)"),
  make_option(c("-s", "--start"), type="numeric", default=0, help="starting value for evolution (root value)"),
  make_option(c("-v", "--verbose"), action="store_false", default=TRUE, help="Print extra output"),
  make_option(c("-z", "--Description"), action="store_false", default=TRUE, help="Script description: simulate continuous trait evolution across a tree")
)
# get command line options, if help option encountered print help and exit, # otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))


### I/O error check
if(is.null(opt$tree)){ stop(" ERROR: provide a tree file (-t)")}

### Functions
rTraitContW = function(phy, weight=1, sigma=0.1, ancestor = FALSE, root.value = 0, ...){
# implementing method from: 
  # author = {Münkemüller, Tamara and Lavergne, Sébastien and Bzeznik, Bruno and Dray, Stéphane and Jombart, Thibaut and Schiffers, Katja and Thuiller, Wilfried},
  # title = {How to measure and test phylogenetic signal},
  # journal = {Methods in Ecology and Evolution}  
# returning values as a dataframe

  # BM = brownian motion
  # weighting formula: w * traitBM + (1-w)*traitRand
  # followed by z-standardizing
  
  require(ape)
  
  ## traitBM
  traitBM <- rTraitCont(phy=phy, ancestor=ancestor, sigma=sigma, root.value=root.value) 
  
  ## traitRand
  traitRand <- sample(traitBM)
  
  ## weighting (1 = all traitBM; 0 = all traitRand)
  x <- weight * traitBM + (1-weight) * traitRand
  
  ## z-standardizing
  x <- scale(x, center=TRUE, scale=TRUE)
  
  return( as.data.frame(x) )
}



### data processing
# loading tree #
if(opt$format == "newick"){ tree <- read.tree(opt$tree) } else
  if(opt$format == "nexus"){ tree <- read.nexus(opt$tree) }

# trait evolution
#data(bird.orders)
res <- rTraitContW(phy=tree, weight=opt$weight, root.value=opt$start)

print(res)

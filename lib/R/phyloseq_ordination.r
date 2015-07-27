#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: phyloseq_ordination.r [options] <phyloseq> <outFile>

options:
  <phyloseq>     Phyloseq object file.
  <outFile>      Name of output file.
                 [Default: ordination.pdf]
  -m=<m>         Ordination method.
                 [Default: NMDS]
  -d=<d>         Distance method.
                 [Default: bray]
  --width=<w>    Plot width.
                 [Default: 8]
  --height=<h>   Plot height.
                 [Default: 6]
  -h             Help

Description:
  Create an ordination from the phyloseq object.
' -> doc

opts = docopt(doc)


# packages
pkgs <- c('phyloseq', 'ggplot2')
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


## making deseq2 object
ord = ordinate(physeq, method=opts[['-m']], distance=opts[['-d']])
p = plot_ordination(physeq, ord, justDF=T)
p$library = as.character(p$library)

ggplot(p, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill=library, size = BD_mid), pch=21) +
    scale_size(range=c(2,8)) +
      theme_bw() +
        theme(
          text = element_text(size=16)
        )

plot.width = as.numeric(opts[['--width']])
plot.height = as.numeric(opts[['--height']])
ggsave(opts[['<outFile>']], width=plot.width, height=plot.height)

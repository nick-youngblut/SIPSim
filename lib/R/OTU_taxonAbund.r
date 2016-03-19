#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: OTU_taxonAbund.r [options] <OTU> 

options:
  <OTU>          OTU table file.
  -o=<o>         Base name of output file.
                 [Default: taxonAbund]
  -t=<t>         File with a list of taxa to select
                 (1 taxon per line). 
  -r=<r>         Max abundance rank to plot (0 = all).
                 [Default: 0]
  -l             Include a legend in the plot.
  --width=<w>    Plot width.
                 [Default: 8]
  --height=<h>   Plot height.
                 [Default: 6]
  -h             Help

Description:
  Plot the abundances of taxa in the OTU table.

  Filtering:
    Taxa filtered by rank (-r) prior to the list (-l).
' -> doc

opts = docopt(doc)


# packages
pkgs <- c('dplyr', 'tidyr', 'ggplot2')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}

# main
## import
tbl = read.delim(opts[['<OTU>']], sep='\t')

## BD min/max/mid
calc_BD_mid = function(BD_min, BD_max){
  if(is.infinite(BD_min)){
    return(BD_max)
  } else
    if(is.infinite(BD_max)){
      return(BD_min)
    } else {
      return((BD_max - BD_min) / 2 + BD_min)
    }
  stop('LOGIC error')
}
#tbl$BD_mid = mapply(calc_BD_mid,tbl$BD_min, tbl$BD_max) 

calc_BD_width = function(BD_min, BD_max){
  if(is.infinite(BD_min)){
    return(0.001)
  } else
    if(is.infinite(BD_max)){
      return(0.001)
    } else {
      return(BD_max - BD_min)
    }
  stop('LOGIC error')
}
tbl$BD_width = mapply(calc_BD_width,tbl$BD_min, tbl$BD_max)


## getting relative abundances (counts)
tbl.s = tbl %>%
  group_by(library, BD_mid) %>%
    summarize(total_sample_count = sum(count))

tbl = inner_join(tbl, tbl.s, c('library'='library','BD_mid'='BD_mid'))
tbl = tbl %>%
  mutate(rel_count = count / total_sample_count)


## filtering / selecting
### selecting abundance ranks of interest
rank.sel = as.numeric(opts[['-r']])
if(rank.sel > 0){
  tbl.s = tbl %>%
    group_by(library, taxon) %>%
      summarize( total_abund = sum(count) ) %>%
        mutate(rank = min_rank(desc(total_abund))) %>%
          filter(rank <= rank.sel)
  tbl = tbl %>%
    filter(taxon %in% tbl.s$taxon)
}
### selecting from taxa from
if(! is.null(opts[['-t']])){
  taxa.to.keep = read.table(opts[['t']], header=FALSE, sep='\t')
  tbl = tbl %>%
    filter(taxon %in% taxa.to.keep[,1])
}

## plotting
BD.GCp0 = 0 * 0.098 + 1.66
BD.GCp50 = 0.5 * 0.098 + 1.66
BD.GCp100 = 1 * 0.098 + 1.66

## making plots
make_frac_plot = function(tbl, BD.GCp50, rel=FALSE, legend=FALSE){
  p = ggplot(tbl, aes(x=BD_mid, fill=taxon))
  if(rel==TRUE){
    p = p + geom_area(aes(y=rel_count), stat='identity', alpha=1, position='stack') +
      labs(y='Relative abundance')
  } else {
    p = p + geom_area(aes(y=count), stat='identity', alpha=0.5, position='dodge') +
      labs(y='Absolute abundance')    
  }  
  p = p + geom_vline(xintercept=c(BD.GCp50), linetype='dashed', alpha=0.5) +
      labs(x='Buoyant density') +
        scale_x_continuous(expand=c(0.01,0)) +
          scale_y_continuous(expand=c(0,0.01)) +
            facet_grid(library ~ .) +
              theme_bw()

  if(legend==TRUE){
    p = p + theme(
      text = element_text(size=16)
    )
  } else {
    p = p + theme(
      text = element_text(size=16),
      legend.position = 'none'
    )
  }
    
  return(p)
}

p.dodge = make_frac_plot(tbl, BD.GCp50, legend=opts[['-l']])
p.fill = make_frac_plot(tbl, BD.GCp50, legend=opts[['-l']], rel=TRUE)


### writing plots
plot.width = as.numeric(opts[['--width']])
plot.height = as.numeric(opts[['--height']])

outFile = paste(c(opts[['-o']], 'abs-abund.pdf'), collapse='_')
ggsave(outFile, plot=p.dodge, width=plot.width, height=plot.height)
message('File written: ', outFile)

outFile = paste(c(opts[['-o']], 'rel-abund.pdf'), collapse='_')
ggsave(outFile, plot=p.fill, width=plot.width, height=plot.height)
message('File written: ', outFile)

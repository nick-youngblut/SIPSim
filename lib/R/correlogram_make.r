#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: correlogram_make.r [options] <data>

options:
  <data>    Either a SIPSim OTU table or a list
            of phyloseq objects. 
  -l        The `data` object is a list of phyloseq objects.
  -h        Help

description:
  Calulcate a correlogram assessing community compostion
  autocorrelations as a function on difference in
  gradient buoyant density. 

  Input should either be 1) a SIPSim OTU table
  from a single SIPSim simulation 2) a list of phyloseq
  objects (eg., communities from multiple days).

  The output is written to STDOUT.
' -> doc

opts = docopt(doc)


# packages
pkgs <- c('dplyr', 'tidyr', 'phyloseq')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}



#-- functions --#

min_max_BD = function(){
  ## min G+C cutoff
  min_GC = 13.5
  ## max G+C cutoff
  max_GC = 80
  ## max G+C shift
  max_13C_shift_in_BD = 0.036

  
  min_BD = min_GC/100.0 * 0.098 + 1.66
  max_BD = max_GC/100.0 * 0.098 + 1.66
  max_BD = max_BD + max_13C_shift_in_BD
  
#  cat('Min BD:', min_BD, '\n')
#  cat('Max BD:', max_BD, '\n')
  return(c('min_BD' = min_BD, 'max_BD' = max_BD))
}


load_simulated = function(filename){
  sim = read.delim(filename, sep='\t') %>%
    select(library, fraction, taxon, BD_mid, rel_abund) %>%
      rename('OTU' = taxon,
             'Buoyant_density' = BD_mid,
             'abundance' = rel_abund,
             'sample' = fraction)
  return(sim)
}

otu2df = function(x){
  x %>% otu_table %>% as.data.frame
}
sample2df = function(x){
  x %>% sample_data %>% as.data.frame
}

load_emperical = function(filename){
  # import object
  emp = readRDS(filename)

  # getting all sample data
  emp_sample_data = do.call(rbind, lapply(emp, sample2df))

  # convert to dataframes
  tmp = lapply(emp, otu2df)

  samps = names(tmp)
  emp = tmp[[samps[1]]]
  emp$OTU = rownames(emp)
  for (x in samps[2:length(samps)]){
        y = tmp[[x]]
            y$OTU = rownames(y)
            emp = left_join(emp, y, c('OTU' = 'OTU'))
      }
  tmp = NULL

  # format dataframe
  emp[is.na(emp)] = 0
  emp = emp %>%
        gather(sample, abundance, starts_with('12C-Con')) %>%
          mutate(sample = sample %>% as.vector)
  emp_sample_data = emp_sample_data %>%
    mutate(X.Sample = X.Sample %>% as.character %>% as.vector)
  
  emp = inner_join(emp, emp_sample_data, c('sample' = 'X.Sample')) %>%
    mutate(Day = Day %>% as.character) %>%
      group_by(sample) %>%
        ungroup() %>%
          select(Day, sample, OTU, Buoyant_density, abundance) %>%
            rename('library' = Day)   # library by day
  return(emp)
}



# delta BD distance matrix
BD.diffs = function(df){
  BDs = df$Buoyant_density %>% unique
  df.BD = expand.grid(BDs, BDs)
  df.BD$diff = df.BD %>% apply(1, diff) %>% abs %>% as.vector
  df.BD = df.BD %>%
    spread(Var1, diff)
  rownames(df.BD) = df.BD$Var2
  df.BD$Var2 = NULL
  dist.BD = df.BD %>% as.matrix
  dist.BD[upper.tri(dist.BD, diag=TRUE)] = 0
  dist.BD %>% as.dist
}

# community distance matrix
vegdist.by = function(df, ...){
  df.w = df %>%
    select(OTU, sample, abundance) %>%
      spread(OTU, abundance) %>%
        as.data.frame()
  rownames(df.w) = df.w$sample
  df.w$sample = NULL
  
  vegan::vegdist(df.w, ...)
}

# correlogram
m.corr = function(X, D, ...){
  res = list()
  for (i in 1:length(X)){
    tmp = vegan::mantel.correlog(X[[i]], D[[i]], ...)
    tmp = tmp['mantel.res'][['mantel.res']] %>% as.data.frame
    colnames(tmp) = c('class.index', 'n.dist', 'Mantel.corr', 'Pr', 'Pr.corr')
    res[[i]] = tmp
  }
  return(res)
}


# calc.corr
calc.corr = function(df, min_BD, max_BD){
  df.corr = df %>%
    filter(Buoyant_density >= min_BD, Buoyant_density <= max_BD) %>%
      group_by(library) %>%
        nest() %>%
          mutate(dist.bray = lapply(data, vegdist.by),
                 dist.BD = lapply(data, BD.diffs),
                 mantel.corr = m.corr(dist.bray, dist.BD, n.class=24)) %>%
                   select(library, mantel.corr) %>%
                     unnest(mantel.corr %>% purrr::map(function(x) x))
  return(df.corr)
}



#--  main --#
BD = min_max_BD()

# loading data
if(opts[['-l']] == TRUE){
  df = load_emperical(opts[['<data>']])
} else {
  df = load_simulated(opts[['<data>']])
}

# calculating correlograms
df.corr = calc.corr(df, BD[1], BD[2])

# writing
write.table(df.corr, stdout(), quote=FALSE, row.names=FALSE, sep='\t')


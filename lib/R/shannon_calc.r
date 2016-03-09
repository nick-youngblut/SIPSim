#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: shannon_calc.r [options] <data> 

options:
  <data>           Either a SIPSim OTU table or a list
                   of phyloseq objects.
  -l               The `data` & `data_preFrac` object
                   are lists of phyloseq objects.
  -h               Help

description:
  Calculate the Shannon index for each gradient fraction community.

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


emp2df = function(x){
  # convert to dataframes
  tmp = lapply(x, otu2df)

  samps = names(tmp)
  emp = tmp[[samps[1]]]
  emp$OTU = rownames(emp)
  for (x in samps[2:length(samps)]){
    y = tmp[[x]]
    y$OTU = rownames(y)
    emp = left_join(emp, y, c('OTU' = 'OTU'))
  }
  tmp = NULL
  emp[is.na(emp)] = 0
  
  return(emp)
}


load_emperical = function(filename){
  # import object
  emp = readRDS(filename)

  # getting all sample data
  emp_sample_data = do.call(rbind, lapply(emp, sample2df))

  # converting to dataframe
  emp = emp2df(emp)
  
  
  # format dataframe
  emp = emp %>%
    gather(sample, abundance, starts_with('12C-Con'))
  emp = inner_join(emp, emp_sample_data, c('sample' = 'X.Sample')) %>%
    mutate(Day = Day %>% as.character) %>%
      group_by(sample) %>%
        ungroup() %>%
          select(Day, sample, OTU, Buoyant_density, abundance) %>%
            rename('library' = Day)   # library by day
  
  return(emp)
}





#--  main --#
BD = min_max_BD()

# loading data
if(opts[['-l']] == TRUE){
  df = load_emperical(opts[['<data>']])
} else {
  df = load_simulated(opts[['<data>']])
}


shannon.long = function(df, abundance_col, ...){
  # calculating shannon diversity index from a 'long' formated table
  ## community_col = name of column defining communities
  ## abundance_col = name of column defining taxon abundances
  df = df %>% as.data.frame
  cmd = paste0(abundance_col, '/sum(', abundance_col, ')')
  df.s = df %>%
    group_by_(...) %>%
      mutate_(REL_abundance = cmd) %>%
        mutate(pi__ln_pi = REL_abundance * log(REL_abundance),
               shannon = -sum(pi__ln_pi, na.rm=TRUE)) %>%
                 ungroup() %>%
                   dplyr::select(-REL_abundance, -pi__ln_pi) %>%
                     distinct_(...)
  return(df.s)
}

# calculating shannon index
df.shan = shannon.long(df, 'abundance', 'library', 'sample') %>%
  filter(Buoyant_density >= BD[1], Buoyant_density <= BD[2]) %>%
  select(-abundance)


# writing
write.table(df.shan, stdout(), quote=FALSE, row.names=FALSE, sep='\t')


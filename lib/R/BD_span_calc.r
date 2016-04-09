#!/usr/bin/Rscript

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'usage: BD_span_calc.r [options] <data> <data_preFrac>

options:
  <data>           Either a SIPSim OTU table or a list
                   of phyloseq objects.
  <data_preFrac>   The pre-fraction data associated with
                   the `data` object.
  -l               The `data` & `data_preFrac` object
                   are lists of phyloseq objects.
  -h               Help

description:
  Calulcate the span of taxa across the buoyant density
  range of the gradient.

  Input should either be 1) a SIPSim OTU table & community file
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


load_simulated_preFrac = function(filename){
  sim.preFrac = read.delim(filename, sep='\t') %>%
    rename('preFrac_abund' = rel_abund_perc,
           'OTU' = taxon_name) %>%
             mutate(preFrac_abund = preFrac_abund / 100) %>%
               dplyr::select(-rank)
  return(sim.preFrac)
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


load_emperical_preFrac = function(filename){
  emp.preFrac = readRDS(filename)

  # getting all sample data
  emp.pf_sample_data = do.call(rbind, lapply(emp.preFrac, sample2df))

  # converting to dataframe
  emp.preFrac = emp2df(emp.preFrac)

  # formatting
  emp.preFrac = emp.preFrac %>%
    gather(sample, abundance, starts_with('12C-Con')) %>%
      mutate(Day = gsub('.+\\.D([0-9]+)\\.R.+', '\\1', sample)) %>%
        rename('preFrac_abund' = abundance,
               'preFrac_sample' = sample,
               'library' = Day)

  return(emp.preFrac)
}



BD_span = function(df){
  df %>%
    filter(abundance > 0) %>%
      group_by(library) %>%
        mutate(max_BD_range = max(Buoyant_density) - min(Buoyant_density)) %>%
          ungroup() %>%
            group_by(OTU, library) %>%
              summarize(mean_preFrac_abund = mean(preFrac_abund),
                        min_BD = min(Buoyant_density),
                        max_BD = max(Buoyant_density),
                        BD_range = max_BD - min_BD,
                        BD_range_perc = BD_range / first(max_BD_range) * 100) %>%
                          ungroup()
}


#--  main --#
n.tile = opts[['-b']] %>% as.numeric
BD = min_max_BD()

# loading data
if(opts[['-l']] == TRUE){
  df = load_emperical(opts[['<data>']])
  df.preFrac = load_emperical_preFrac(opts[['<data_preFrac>']])
} else {
  df = load_simulated(opts[['<data>']])
  df.preFrac = load_simulated_preFrac(opts[['<data_preFrac>']])
}


# joining
df.j = inner_join(df, df.preFrac, c('OTU' = 'OTU', 'library' = 'library'))


# calculating BD span
df.j.BDspan = BD_span(df.j)


# writing
write.table(df.j.BDspan, stdout(), quote=FALSE, row.names=FALSE, sep='\t')


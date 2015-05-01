#!/usr/bin/env python

#--- Option parsing ---#
"""
isoIncorp: set taxon isotope incorporation distributions

Usage:
  isoIncorp [options] <BD_KDE> <config_file>
  isoIncorp -h | --help
  isoIncorp --version

Options:
  <BD_KDE>           Buoyant Density KDE object file ('-' if from STDIN).
  <config_file>      Config file setting inter-population-level incorporation
                     distribution (see Description).
  --comm=<cf>        Community file used set abundance-weighting of isotope.
  -n=<n>             Number of Monte Carlo replicates to estimate
                     G+C error due to diffusion. 
                     [default: 100000]
  --isotope=<is>     Isotope incorporated by taxa (13C or 15N).
                     [default: 13C]
  --phylo=<phy>      Newick phylogeny of taxa used for brownian motion evolution 
                     of distribution parameters.
  --bw=<bw>          The bandwidth scalar or function passed to
                     scipy.stats.gaussian_kde().
  --np=<np>          Number of parallel processes.
                     [default: 1]
  --cs=<cs>          Chunksize for each process (number of taxa).
                     [default: 1]
  -h --help          Show this screen.
  --version          Show version.
  --debug            Debug mode

Description:
  For each population (taxon), there is a distribution of how much
  isotope is incorporated by each member of the population (intra-population
  isotope incorporation).
  These intra-population distributions are, in reality, likely to
  vary *among* populations (inter-population variation). This command is used
  to set the intra-population distributions of isotope abundance
  (e.g., normal or uniform) and how the parameters of those distributions
  (e.g., the mean or sd) vary *among* taxa (inter-population).

  In other words, for each population (taxon), you are setting the distribution of isotope
  incorporation for that population, and these population-level distributions 
  can vary among populations (taxa).

  For example, population-level incorporation is drawn from a normal distribution
  for each taxon, but the mean (mu) of each population-level distribution is
  determined by a uniform distribution that varies from 0 to 100. Thus, while
  all populations have a normal distribution of incorporation, some will have
  much more incorporation (near 100%) than others (close to 0%).

  Distributions can be standard distributions (eg., normal or uniform) or a
  mixture model, which is a combination of weighted distributions.
  These distribution parameters can be defined with a config file (see below).

  Inter-population variation in isotope incorporation can alternatively be 
  modeled as Brownian motion evolution. In this case, each intra-population
  parameter (eg., the 'mu' parameter for a normal distribution) is evolved
  across the user-provided phylogeny as Brownian motion evolution.
  WHY USE THIS METHOD? Brownian motion evoution will generally give highly 
  related taxa more similar values than more distantly related taxa. This may
  be more biologically realistic than simulating inter-population variation
  as a normal or other standard distribution. 
  As a null model of Brownian motion evolution, random evolution can be set
  instead (no correlation between taxon relatedness and evolved character
  value similarity). These 2 models of evolution can be combined with different
  contributes ('ratio' parameter in the config file).


  Supported intra-population distributions:
    'normal' = np.random.normal
    'uniform' = np.random.uniform
    'BM' = browian motion evolution (must provide --phylo)

  OUTPUT:
    An updated BD KDE object file is written to STDOUT.

  NOTES:
    * Multiple standard distriutions can be provided to create a mixture model.

  CONFIG FILE:
    File format:  http://www.voidspace.org.uk/python/configobj.html#config-files
    * NOTE: the goal is to set INTRA-population parameters, which is why 
      intra-population is higher in the hierarchy than inter-population
  ----------------- START ------------------
  [library1]
    max_perc_taxa_incorp = 100   # <- max percent of taxa with any isotope incorporation

    [[intraPopDist1]]  # <- the intra-population isotope distrubution
      distribution = normal

      [[[loc]]]     # <- the inter-pop variation in 'loc' param for intra-pop distribution
        
        [[[[interPopDist 1]]]]   # <- you can have >1 standard distribution 
        distribution = normal
        loc = 90
        scale = 2

      [[[scale]]]
        
        [[[[interPopDist 1]]]]
        distribution = normal
        loc = 10
        scale = 2

  -----------------  END   ------------------

     * Multiple intraPopDist or interPopDist:
       * If multiple distributions provided, then each will be incorporated
         into a mixture model.
       * A 'weight' parameter should be provided 
          * example: weight = 0.5    # 50% weighting 

     * Using Brownian Motion evolution (must provide --phylo):
       * params:
         * ratio = ratio of brownian motion vs random sampling used 
         * sigma = standard deviation for evolving values
                   (drawing from a normal distribution)
"""

# import
## batteries
from docopt import docopt
import sys,os
from functools import partial
import re
import cPickle as pickle
from random import shuffle

## 3rd party
import pandas as pd
import parmap
import scipy.stats as stats

## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

##import IsoIncorp
from SIPSim import CommTable
import IsoIncorp
import Utils
from SIPSimCython import add_incorp


# functions
def _make_kde(taxon_name, x, libID, config, taxa_incorp_list, 
             isotope='13C', n=10000, bw_method=None): 
    """Making new KDE of BD value distribution which includes
    BD shift due to isotope incorporation. 
    Args:
    taxon_name -- str; name of taxon
    x -- dict; keys: kde, [abundances]
    libID -- str; library ID
    config -- config object
    taxa_incorp_list -- iterable; taxa that can incorporate isotope
    isotope -- str; isotope that is incorporated
    n -- number of Monte Carlo samples to use for estimating
         BD+isotope_BD distribution
    bw_method -- bandwidth scalar or function passed to
                 scipy.stats.gaussian_kde().
    Return:
    (taxon_name, KDE*)
       * Note: KDE object may be None    
    """

    # status
    sys.stderr.write('Processing: {}\n'.format(taxon_name))

    # unpack
    n = int(n)
    kde = x['kde']
    if kde is None:
        return (taxon_name, None)

    # can taxon incorporate any isotope?
    if taxon_name not in taxa_incorp_list:
        return (taxon_name, None)

    # taxon abundance for library
    try:
        taxon_abund = x['abundances'][libID]
    except KeyError:
        tmp = re.sub('[Ll]ib(rary)* *','', libID)
        try:
            taxon_abund = x['abundances'][tmp]
        except KeyError:
            taxon_abund = None

    # max incorporation for isotope
    maxIsotopeBD = IsoIncorp.isotopeMaxBD(isotope)

    # making a mixture model object lib:taxon
    ## mix_model => distribution of % incorporation for taxon population
    mix_model = IsoIncorp.make_incorp_model(taxon_name, libID, config)

    # making KDE of BD + BD_isotope_incorp
    kdeBD = stats.gaussian_kde( 
        add_incorp(kde.resample(n)[0], 
                   mix_model, 
                   maxIsotopeBD),
        bw_method=bw_method)

    return (taxon_name, kdeBD)


def _add_comm_to_kde(KDE_BD, comm):
    """Adding comm data for each taxon to each KDE.
    Args:
    Return:
    KDE_BD -- {taxon_name:{kde|abundances}}
    """
    libIDs = comm.get_unique_libIDs()
    for taxon_name,kde in KDE_BD.items():
        d = {'kde':kde, 'abundances':{}}
        for libID in libIDs:
            abund = comm.get_taxonAbund(taxon_name, libID=libID)
            d['abundances'][libID] = abund[0]
        KDE_BD[taxon_name] = d    


def _taxon_incorp_list(libID, config, KDE_BD):
    # perc incorp from config
    try:
        max_perc_taxa_incorp = config.get_max_perc_taxa_incorp(libID) 
    except KeyError:
        max_perc_taxa_incorp = 100.0
    max_perc_taxa_incorp /= 100.0

    # randomized list of taxa
    taxon_names = KDE_BD.keys()
    shuffle(taxon_names)

    # set of taxa that incorporate any isotope
    n_incorp = int(round(len(taxon_names) * max_perc_taxa_incorp, 0))
    return taxon_names[:n_incorp]



def main(args):
    # loading input
    ## kde_bd
    KDE_BD = Utils.load_kde(args['<BD_KDE>'])
    ## comm (optional)
    if args['--comm'] is not None:
        comm = CommTable.from_csv(args['--comm'], sep='\t')
    else:
        comm = None


    # combining kde and comm
    _add_comm_to_kde(KDE_BD, comm)
    
    # loading the config file
    config = IsoIncorp.Config.load_config(args['<config_file>'],
                                          phylo=args['--phylo'])
     

    # creating kde of BD distributions with BD shift from isotope
    KDE_BD_iso = dict()
    for libID in config.keys():        
        # making a list of taxa that can incorporate 
        taxa_incorp_list = _taxon_incorp_list(libID, config, KDE_BD)

        # setting params for parallelized function
        pfunc = partial(_make_kde, 
                        config = config,
                        libID = libID,
                        n = args['-n'],
                        taxa_incorp_list = taxa_incorp_list,
                        isotope = args['--isotope'],
                        bw_method = args['--bw'])                    

        # parallel by taxon
        tmp = parmap.starmap(pfunc, KDE_BD.items(),
                             processes = int(args['--np']),
                             chunksize = int(args['--cs']),
                             parallel = not args['--debug'])

        KDE_BD_iso[libID] = {taxon:KDE for taxon,KDE in tmp}
    

    # writing pickled BD-KDE with BD shift from isotope incorp
    pickle.dump(KDE_BD_iso, sys.stdout)
        
    
    
# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')    
    main(args)



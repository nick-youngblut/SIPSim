#!/usr/bin/env python

#--- Option parsing ---#
"""
isoIncorp: set taxon isotope incorporation distributions

Usage:
  isoIncorp [options] <comm_file> <config_file>
  isoIncorp -h | --help
  isoIncorp --version

Options:
  <comm_file>          Output from `SIPSim gradientComm`. ('-' if from STDIN)
  <config_file>        Config file setting inter-population-level incorporation
                       distribution (see Description).
  --phylo=<phy>        Newick phylogeny of taxa used for brownian motion evolution 
                       of distribution parameters (NOT YET IMPLEMENTED).
  -h --help            Show this screen.
  --version            Show version.
  --debug              Debug mode

Description:
  For each population (taxon), there is a distribution of how much
  isotope is incorporated by each member of the population.
  These intra-population distributions are, in reality, likely to
  vary among populations (inter-population). This command is used
  to set the intra-population distributions of isotope abundance
  (e.g., normal or uniform) and how the parameters of those distributions
  (e.g., the mean or sd) vary among taxa.

  This inter-population variation can be determined either by:
    * Simulating brownian motion evolution of the intra-pop distribution parameters.
      * For example: the intra-pop mean is 'evolved' across the phylogeny, resulting
        in generally more similar intra-pop means among more closely related taxa.
      * The params to be 'evolved' are provided via '--evoStart'
      * The degree of Browian motion vs random selection of values is controlled 
        by '--evoWeight'
    * Selecting the inter-pop isotope incorporation distribution params from user-defined
      inter-pop distributions describing how those params vary among populations.
      * This is set with a configure file (see 'config file' below)
      * A gaussian mixture model can be used to model populations that are
        subdivided in how members incorporate isotope (e.g., high mean incorp in Group 1,
        and low mean incorp in Group 2)

  Supported intra-population incorp distributions (--popIncorp):
    'normal' = np.random.normal
    'uniform' = np.random.uniform

  Config file:
    File format:  http://www.voidspace.org.uk/python/configobj.html#config-files
    * NOTE: the goal is to set INTRA-population parameters, which is why 
      intra-population is higher in the hierarchy than inter-population
  ----------------- START ------------------
  [library1]
    max_percent_incorp = 100   # <- max percent of taxa with any isotope incorporation

    [[intraPopDist1]]  # <- the intra-population isotope distrubution
      distribution = normal

      [[[loc]]]     # <- the inter-pop variation in 'loc' param for intra-pop distribution
      distribution = normal
      loc = 90
      scale = 2

      [[[scale]]]
      distribution = normal
      loc = 10
      scale = 2
  -----------------  END   ------------------

     * Multiple intraPopDist or interPopDist:
       * If multiple distributions provided, then each will be incorporated
         into a mixture model.
       * A 'weight' parameter should be provided.

     * Using Brownian Motion evolution (must provide --phylo):
       * params:
         * root = value assigned for root which will then be evolved across the tree
         * ratio = ratio of brownian motion vs random sampling used 
         * sigma = standard deviation for evolving values (drawing from a normal distribution)

"""

# import
## batteries
from docopt import docopt
import sys,os

## 3rd party
import pandas as pd

## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

#import IsoIncorp
from SIPSim import CommTable
import IsoIncorp


# functions
def main(Uargs):
    """
    Uargs -- dict of user-provided args
    """
    # loading community file as dataframe
    #comm = pd.read_csv(Uargs['<comm_file>'], sep='\t')

    # loading community file
    comm = CommTable.from_csv(Uargs['<comm_file>'], sep='\t')

    # loading the config file
    config = IsoIncorp.Config.load_config(Uargs['<config_file>'],
                                          phylo=Uargs['--phylo'])
    
    distsTbl = IsoIncorp.populationDistributions(config, comm)

    # writing
    distsTbl.to_csv(sys.stdout, sep='\t', index=None)
    
    
    # if --config
    #! set inter-pop distributions based on config {param:numpy_dist}
    ## loading config file and set distributions
    ## foreach sample:
    ### randomly select taxon (up to % of taxa needed):
    #### draw params from inter-pop distributions
    #### write: sample,taxon,distribution,dist_params...
    ### other taxa:
    #### params = uniform dist, min-max=0
    
    # elif --phylo:
    #! calling rTraitPhylo.r
    ## foreach sample:
    ### foreach intra-pop param:
    #### call rTraitPhylo.r; save values
    ### foreach taxon:
    #### apply values to indiv. pop dist
    #### write: sample,taxon,distribution,dist_params...
    ### other taxa:
    #### params = uniform dist, min-max=0
    
    
    
# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')    
    main(Uargs)



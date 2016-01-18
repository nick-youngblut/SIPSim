#!/usr/bin/env python

#--- Option parsing ---#
"""
isotope_incorp: set taxon isotope incorporation distributions

Usage:
  isotope_incorp [options] <BD_KDE> <config_file>
  isotope_incorp -h | --help
  isotope_incorp --version

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

  In other words, for each population (taxon), you are setting the distribution
  of isotope incorporation for that population, and these population-level
  distributions can vary among populations (taxa).

  For example, population-level incorporation is drawn from a normal
  distribution for each taxon, but the mean (mu) of each population-level
  distribution is determined by a uniform distribution that varies from 0 to
  100. Thus, while all populations have a normal distribution of incorporation, 
  some will have much more incorporation (near 100%) than others (close to 0%).

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

  Output
  ------
  An updated BD KDE object file is written to STDOUT.

  Notes
  -----
  * Multiple standard distriutions can be provided to create a mixture model.

  Config file
  -----------
    File format:  http://www.voidspace.org.uk/python/configobj.html#config-files
    * NOTE: the goal is to set INTRA-population parameters, which is why 
      intra-population is higher in the hierarchy than inter-population
    ------ START ------

  # top level is 'library' (gradient); this must be an integer 
  [1]
    max_perc_taxa_incorp = 100  #<- max percent of taxa with any isotope incorp

    [[intraPopDist1]]  # <- intra-population isotope distrubution
      distribution = normal

      [[[loc]]]   # <- inter-pop variation in 'loc' param for intra-pop dist.
        
        [[[[interPopDist 1]]]]   # <- you can have >1 standard distribution 
        distribution = normal
        loc = 90
        scale = 2

      [[[scale]]]
        
        [[[[interPopDist 1]]]]
        distribution = normal
        loc = 10
        scale = 2

    ------  END  ------

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

## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

##import IsoIncorp
import IsoIncorp

    
# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')    
    IsoIncorp.main(args)



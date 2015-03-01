#!/usr/bin/env python

#--- Option parsing ---#
"""
gradientComms: simulate communities in the samples used for SIP

Usage:
  gradientComms [options] <taxon_list>
  gradientComms -h | --help
  gradientComms --version

Options:
  <taxon_list>        A file listing taxon names ('-' if from STDIN). [default: -]
  --shared_perc=<sp>  The percent of taxa shared in each community.
                      Percent set by the community with the smallest richness.
                      Example: if smallest community is 10 taxa,
                               a shared percent of 20% = 2 shared taxa.
                      The total taxon pool must be large enough to accommodate all un-shared taxa.
                      [default: 100]
  --richness=<r>      The number of taxa in each library.
                      Values of 0-1 will be interpreted as a fraction of the total taxa pool.
                      [default: 1]
  --abund_dist=<a>    The statistical distribution used for selecting relative abundances.
                      (see numpy.random for a list of distributions).
                      [default: exponential]
  --abund_dist_p=<p>  Abundance distribution parameters.
                      (see numpy.random for distribution params).
                      [default: scale:1]
  --perm_perc=<pp>    How much to vary the rank-abundances between communities.
                      [default: 0]
  --n_comm=<nc>       Number of communities to simulate.
                      [default: 1]
  --config=<c>        Config file for setting community-specific parameters (& global params).
                      Community-specific parameters can include:
                        * richness
                        * abund_dist
  --debug             Debug mode
  --version           Show version.
  -h --help           Show this screen.

Description:
  Simulating the alpha and and beta diversities of >=1 community.

  Output:
    A tab-delimited table of taxon abundances for each library is written to STDOUT.

"""

# import
## batteries
from docopt import docopt
import sys,os


## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from SimComms import SimComms


# functions
def call_grinder(genomeListFile, filePath=None, profileFile=None, exe='grinder'):
    """System call of grinder
    Args:
    genomeListFile -- genome list file name
    filePath -- full path to genome files (if needed)
    exe -- name of grinder executable
    profileFile -- profile file to provide to grinder
    """

    cmd = '{exe} -rf {glf} -fp {fp}'.format(exe=exe, glf=genomeListFile, fp=filePath)
    
    if profileFile is not None:
        cmd += ' -pf {}'.format(profileFile)
        
#    print subprocess.check_output([cmd], shell=True)
    


def main(uargs):

    # init
    SC = SimComms(taxon_list = uargs['<taxon_list>'],
                  perm_perc = uargs['--perm_perc'],
                  shared_perc = uargs['--shared_perc'],
                  richness = uargs['--richness'],
                  abund_dist = uargs['--abund_dist'],
                  abund_dist_params = uargs['--abund_dist_p'],
                  n_comm = uargs['--n_comm'],
                  config = uargs['--config'])

    # making communities
    for i in xrange(SC.n_comm):
        SC.make_comm(str(i+1))

    
if __name__ == '__main__':
    uargs = docopt(__doc__, version='0.1')
    main(uargs)

    
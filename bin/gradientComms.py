#!/usr/bin/env python

#--- Option parsing ---#
"""
gradientComms: simulate communities in the samples used for SIP

Usage:
  gradientComms [options] <genomeList>
  gradientComms -h | --help
  gradientComms --version

Options:
  <genomeList>        A file listing: taxonName<tab>genomeSeqFileName
  --n_comm=<nc>       Number of communities to simulate.
                      [default: 1]
  --richness=<r>      The number of taxa in each library.
                      Values of 0-1 will be interpreted as a fraction of the total taxa pool.
                      [default: 1]
  --abund_dist=<a>    The statistical distribution used for selecting relative abundances.
                      (see numpy.random for a list of distributions).
                      [default: lognormal]
  --abund_dist_p=<p>  Abundance distribution parameters.
                      (see numpy.random for distribution params).
                      [default: mean:10,sigma:2]
  --shared_perc=<sp>  The percent of taxa shared in each community.
                      Percent set by the community with the smallest richness.
                      Example: if smallest community is 10 taxa,
                               a shared percent of 20% = 2 shared taxa.
                      The total taxon pool must be large enough to accommodate all un-shared taxa.
                      [default: 100]
  --perm_perc=<pp>    How much to vary the rank-abundances between communities.
                      [default: 0]
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
def main(uargs):
    # init
    SC = SimComms(taxon_list = uargs['<genomeList>'],
                  perm_perc = uargs['--perm_perc'],
                  shared_perc = uargs['--shared_perc'],
                  richness = uargs['--richness'],
                  abund_dist = uargs['--abund_dist'],
                  abund_dist_params = uargs['--abund_dist_p'],
                  n_comm = uargs['--n_comm'],
                  config = uargs['--config'])

    # making communities
    for comm_id in SC.keys():        
        SC.make_comm(comm_id)

    # permuting based on permute perc
    for comm_id, comm in SC.items():
        SC.permute(comm, uargs['--perm_perc'])

    # writing out table
    SC.write_comm_table()

        
if __name__ == '__main__':
    uargs = docopt(__doc__, version='0.1')
    main(uargs)

    

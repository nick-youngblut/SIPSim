#!/usr/bin/env python

#--- Option parsing ---#
"""
OTU_sim: simulate OTUs for gradient fractions based on simulated
fragment G+C content and isotope incorporation

Usage:
  OTU_table [options] <BD_KDE> <communities> <fractions>
  OTU_table -h | --help
  OTU_table --version

Options:
  <BD_KDE>            KDE object of BD value distributions. 
                      ('-' if from STDIN)
  <communities>       Simulated community abundance table file.
  <fractions>         Simulated gradient fraction file.
  --abs=<aa>          Absolute abundance of all taxa in the community. 
                      [default: 1e5]
  --quiet             Limit STDERR output.
  --version           Show version.
  --debug             Debug mode
  -h --help           Show this screen.


Description:
  Create an OTU table of simulated OTUs for each fraction in >=1 CsCl gradient.

  Basically, the location within the gradient (i.e., buoyant density)
  of each DNA fragment associated with each taxon is determined, and
  then binned into simulated gradient fractions that span certain
  buoyant density ranges.

  The abundance of each OTU in each fraction is based on:
    1) The absolute abundance of the OTU in the pre-gradient community.
    2) The G+C content of each simulated fragment of the taxon.
    3) The fragment length of each simulated fragment (influences diffusion).
"""

# import
## batteries
from docopt import docopt
import os, sys
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import OTU_Table
    

# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')
    OTU_Table.main(args)
    

        

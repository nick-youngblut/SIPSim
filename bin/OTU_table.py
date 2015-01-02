#!/usr/bin/env python

#--- Option parsing ---#
"""
OTU_table: create an OTU table of gradient fractions based on simulated
fragment G+C content and isotope incorporation

Usage:
  OTU_table [options] <fragGC_file> <comm_file> <incorp_file> <frac_file>
  OTU_table -h | --help
  OTU_table --version

Options:
  <fragGC_file>       Name of file produced by fragGC subcommand.
  <comm_file>         Name of file produced by gradientComms subcommand.
  <incorp_file>       Name of file produced by isoIncorp subcommand.
  <frac_file>         Name of file produced by fractions subcommand.
  --abs_abund=<aa>    Absolute abundance of all taxa in the community. [default: 1e6]
  --g_noise=<gn>      scipy distribution function describing gradient 'noise'. [default: cauchy]
  --gn_scale=<np>     Scale parameter for the '--g_noise' distribution. [default: 0.0]
  --gc_range=<gcr>    Min-max possible G+C values post-diffusion or post-noise. [default: 0,100]
  --a_weight=<aw>     Abundance weighting for isotope incorporation.
  --isotope=<is>      Isotope incorporated by taxa (13C or 15N). [default: 13C]
  -h --help           Show this screen.
  --version           Show version.
  --debug             Debug mode

Description:

"""

# import
## batteries
from docopt import docopt
import os, sys
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import SIPSim 
        
# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    SIPSim.main(Uargs)
    

        
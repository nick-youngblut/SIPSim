#!/usr/bin/env python

#--- Option parsing ---#
"""
OTU_sim: simulate OTUs for gradient fractions based on simulated
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
  --isotope=<is>      Isotope incorporated by taxa (13C or 15N). [default: 13C]
  --a_weight=<aw>     Abundance weighting for isotope incorporation (NOT YET IMPLEMENTED).
  --g_noise=<gn>      scipy distribution function describing gradient 'noise'. [default: cauchy]
  --gn_scale=<np>     Scale parameter for the '--g_noise' distribution. [default: 0.0]
  --gc_range=<gcr>    Min-max possible G+C values post-diffusion or post-noise. [default: 0,100]
  --log=<lg>          File name for log of GC and fragment length values ('None' = no file). [default: None]
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

import SIPSim 
        
# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    SIPSim.main(Uargs)
    

        

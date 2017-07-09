#!/usr/bin/env python

"""
qSIP_atomExcess: calculate isotope enrichment from qSIP data

Usage:
  qSIP_atomExcess [options] <OTU_table> <exp_design>
  qSIP_atomExcess -h | --help
  qSIP_atomExcess --version

Options:
  <OTU_table>     OTU table file 
                  (must contain an 'abs_abund' column).
  <exp_design>    Experimental design table. (See Description)
  -i=<i>          Isotope (13C or 18O).
                  [Default: 13C]
  -n=<n>          Number of bootstrap replicates to calculate CIs.
                  [Default: 1000]
  -a=<a>          Alpha for confidence interval.
                  [Default: 0.1] 
  --np=<np>       Number of processors.
                  [Default: 1]
  --byBoot        Parallelization by bootstrap replicate instead of taxon.
                  (useful if running many bootstrap reps on few taxa)
  --version       Show version.
  --debug         Debug mode (no multiprocessing).
  -h --help       Show this screen.


Description:

  'exp_design' input file
  -----------------------
  2-column table: <library><tab><control|treatment>

  __Example-start__
  1<tab>control
  2<tab>treatment
  3<tab>control
  4<tab>treatment
  --Example-end--

References:
  Hungate BA, Mau RL, Schwartz E, Caporaso JG, Dijkstra P, Gestel N van, et
  al. (2015). Quantitative Microbial Ecology Through Stable Isotope Probing.
  Appl Environ Microbiol AEM.02280-15.
"""

# import
## batteries
from docopt import docopt
import sys
import os
from functools import partial
## 3rd party
import numpy as np
import pandas as pd
## application libraries
#scriptDir = os.path.dirname(__file__)
#libDir = os.path.join(scriptDir, '../lib/')
#sys.path.append(libDir)

from SIPSim.OTU_Table import OTU_table
from SIPSim import QSIP_atomExcess
    

def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
        
    otu = QSIP_atomExcess.qSIP_atomExcess(args)
    otu.to_csv(sys.stdout, sep='\t', index=False)
     

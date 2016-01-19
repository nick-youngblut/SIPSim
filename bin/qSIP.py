#!/usr/bin/env python

#--- Option parsing ---#
"""
qSIP: simulate quantitative SIP data

Usage:
  qSIP [options] <OTU_table> <OTU_subsample_table>
  qSIP -h | --help
  qSIP --version

Options:
  <OTU_table>             OTU table file ('true abundances').
  <OTU_subsample_table>   OTU table file (post-sequencing abundances;
                                          relative abundances).
  --reps=<rp>             Number of qPCR replicates.
                          [Default: 3]
  -f=<f>                  Formula describing qPCR value error.
                          Relating mean values to variance (var ~ mean).
                          [Default: 5889 + 1*x + 0.714*x**2]
  --version               Show version.
  --debug                 Debug mode.
  -h --help               Show this screen.


Description:
  Simulate quantitative stable isotope probing (qSIP) data (Hungate et al.,)
  qPCR copy number values are derived from the absolute ('true') taxon
  abundances in the OTU_table file. Error is added to the qPCR values
  according to the formula describing how variance relates to the mean 
  (var ~ mean). The default function was derived from the qPCR data in 
  Hungate et al., 2015. 

  Edited OTU table output
  -----------------------
  The 'prop_abs_abund' column contains the proportional absolute abundances
  values for each taxon (prop_abs_abund = rel_abund * qPCR_total_copy_number).


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
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import QSIP
import Utils
    

# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    #Uargs['-p'] = Utils.parseKeyValueString(Uargs['-p'])

    otu = QSIP.qSIP(Uargs)
    otu.to_csv(sys.stdout, sep='\t', index=False)

    

        

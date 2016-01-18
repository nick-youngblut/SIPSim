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
  -r=<r>                  Binomial distribution dispersion parameter.
                          [Default: 100]
  --version               Show version.
  --debug                 Debug mode.
  -h --help               Show this screen.


Description:
  Simulate quantitative stable isotope probing (qSIP) data (Hungate et al.,)
  qPCR copy number values are derived from the absolute ('true') taxon
  abundances in the OTU_table file. Error is added to the qPCR values
  by scaling variance with the mean as from a binomial distribution:
  (m_err = m + m**/r), where 'm' = mean, and 'r' = the dispersion parameter
  ('-r'). The lower the '-r' value, the more heteroskedasticity.

  qPCR value table (output)
  -------------------------
  The 'total_count' column contains the total counts in the <otu_table> file.
  The 'total_count_qPCR" column contains the total counts with simulated
  qPCR error. 

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

    

# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    otu = QSIP.qSIP(Uargs)
    otu.to_csv(sys.stdout, sep='\t', index=False)

    

        

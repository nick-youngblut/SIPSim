#!/usr/bin/env python

"""
deltaBD: simulate quantitative SIP data

Usage:
  deltaBD [options] <OTU_table> <exp_design> 
  deltaBD -h | --help
  deltaBD --version

Options:
  <OTU_table>    OTU table file ('true abundances').
  <exp_design>   Experimental design table. (See Description)
  -b=<b>         The number of evenly spaced BD bins for linearly
                 interpolating abundances.
                 [Default: 20]
  --version      Show version.
  --debug        Debug mode.
  -h --help      Show this screen.

Description:
  Calculate delta BD as described in Pepe-Ranney and Campell et al., (2015).

  'exp_design' input file
  -----------------------
  2-column table: <library><tab><control|treatment>

  __Example-start__
  1<tab>control
  2<tab>treatment
  3<tab>control
  4<tab>treatment
  --Example-end--

  OUTPUT
  ------
  A table with center of mass (CM) and delta BD values is written to STDOUT.


References:
  Pepe-Ranney C, Campbell AN, Koechli C, Berthrong ST, Buckley DH. (2015).
  Unearthing the microbial ecology of soil carbon cycling with DNA-SIP.
  bioRxiv 022483.

"""

# import
## batteries
from docopt import docopt
import sys
import os
from functools import partial
## 3rd party
## application libraries
from SIPSim import DeltaBD
    


def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
          
    df_deltaBD = DeltaBD.deltaBD(args)
    df_deltaBD.to_csv(sys.stdout, sep='\t', index=False)


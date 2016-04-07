#!/usr/bin/env python

#--- Option parsing ---#
"""
HRSIP: conduct high-resolution SIP method 

Usage:
  HRSIP [options] <phyloseq>
  HRSIP -h | --help
  HRSIP --version

Options:
  <phyloseq>     Phyloseq object. 
  -w=<w>         "heavy" BD window(s). (omma separated list)
                 [Default: 1.71-1.75]
  --cont=<c>     Control libraries. (comma-separated list)
                 [Default: 1]
  --treat=<t>    Treatment libraries. (comma-separated list)
                 [Default: 2]
  --log2=<l>     Log2 fold change cutoff.
                 [Default: 0.25]
  --hypo=<h>     altHypothesis tested by DESeq
                 ("greaterAbs","greater","less")
                 [Default: greater]
  --prefix=<p>   Output file prefix. (Use "None" to use phyyloseq object name).
                 [Default: None]
  --hier         Hierarchical multi-window HR-SIP (HMW-HR-SIP). 
  --padj=<pa>    Adjusted P-value cutoff for calling incorporators.
                 NOTE: only used for hierarchial HR-SIP.
                 [Default: 0.1]
  --version      Show version.
  --debug        Debug mode.
  -h --help      Show this screen.


Description:
  Conduct HR-SIP method for identifying incorporators as detailed in Pepe-Ranney
  & Campbell et al., (in review). 

  `hier`option
  ------------
  This will run hierarchical multi-window HR-SIP (HMW-HR-SIP). This differs
  from multi-window HR-SIP (MW-HR-SIP) in that, for each window, taxa 
  identified as incorporators are filtered out before the next window. This
  reduces the total number of hypotheses that need to be corrected for. 
  The order of windows processed is heaviest to lightest. 

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
## 3rd party
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
import HR_SIP
    

# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    HR_SIP.HR_SIP(Uargs)


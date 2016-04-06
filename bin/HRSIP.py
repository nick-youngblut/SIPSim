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
                 [Default:  None]
  --version      Show version.
  --debug        Debug mode.
  -h --help      Show this screen.


Description:
  Conduct HR-SIP method for identifying incorporators as detailed in Pepe-Ranney
  & Campbell et al., (in review). 

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
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import Utils
    

def parse_BD_windows(x):
    return [y.split('-') for y in x.split(',')]


def HRSIP_by_window(BD_min, BD_max, prefix, Uargs):
    """HRSIP pipeline wrapper
    """
    # trimming physeq to just fractions in BD window
    exe = os.path.join(libDir, 'R', 'phyloseq_edit.r')
    inFile = Uargs['<phyloseq>']
    editFile = prefix + '_{}-{}.physeq'.format(BD_min, BD_max)
    cmd = '{} --BD_min {} --BD_max {} {} > {}'
    cmd = cmd.format(exe, BD_min, BD_max,inFile, editFile)
    Utils.sys_call(cmd)

    # calling DESeq2                                                      
    exe = os.path.join(libDir, 'R', 'phyloseq_DESeq2.r')
    outFile = os.path.splitext(editFile)[0] + '_DESeq2'
    cmd = '{} --log2 {} --hypo {} --cont {} --treat {} {} > {}'
    cmd = cmd.format(exe, Uargs['--log2'], Uargs['--hypo'], Uargs['--cont'],
                     Uargs['--treat'], editFile, outFile)
    Utils.sys_call(cmd)
    
    # returning output file
    return outFile


def combine_DESeq(files, prefix):
    """Combining multiple DESeq2 objects; global p-value adjustment
    """
    exe = os.path.join(libDir, 'R', 'DESeq2_combine.r')
    outFile = prefix + '_DESeq2'
    cmd = '{} {} > {}'.format(exe, ' '.join(files), outFile)
    Utils.sys_call(cmd)
    return outFile


# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')

    # file prefix
    if Uargs['--prefix'].lower() == 'none':
        prefix = os.path.splitext(Uargs['<phyloseq>'])
    else:
        prefix = Uargs['--prefix']

    # parsing BD-windows
    BD_windows = parse_BD_windows(Uargs['-w'])
    
    # DESeq2 on each BD window
    res_files = []
    for BD_min,BD_max in BD_windows:        
        msg = 'HR-SIP on BD window: {}-{}\n'.format(BD_min, BD_max)
        sys.stderr.write(msg)
        resFile = HRSIP_by_window(BD_min, BD_max, prefix, Uargs)
        res_files.append([resFile, BD_min, BD_max])

    # combining DESeq2 files (if multiple)
    if len(res_files) > 1:
        outFile = combine_DESeq([x[0] for x in res_files], prefix)
    else:
        print 'Combining DESeq objects'
        outFile = res_files[0][0]
        outFileNew = prefix + '_DESeq2'
        os.rename(outFile, outFileNew)

    # status
    print 'File written: {}'.format(outFile)

    

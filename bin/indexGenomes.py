#!/usr/bin/env python

#--- Option parsing ---#
"""
indexGenomes: index genomes for in-silico PCR; required for amplicon-fragment simulation

Usage:
  indexGenomes [options] <genomeList>
  indexGenomes -h | --help
  indexGenomes --version

Options:
  <genomeList>    A file listing: taxonName<tab>genomeSeqFileName
  --fp=<fp>       Full path to genomeSeqFiles (if not in genomeList file).
  --K_value=<kv>  Kvalue for indexing. [Default: 9]
  --np=<np>       Number of genomes to process in parallel. [Default: 1]
  --quiet         Limit stderr output. 
  -h --help       Show this screen.
  --version       Show version.
  --debug           Debug mode (no parallel processes)

Description:
  The genomeList file (tab-delim; 'taxon_name<tab>file_name')
  should list all genome sequence files
  (fasta file format; 1 genome per file; genomes can be multi-chromosome/scaffold).

  This script is not needed if shotgun-fragments are to be simulated
  (instead of amplicon-fragments).

"""

# import
## batteries
from docopt import docopt
import sys,os
## 3rd party
## application libraries
scriptDir = os.path.abspath(os.path.dirname(__file__))
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import IndexGenomes
    

# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    IndexGenomes.main(Uargs)
    

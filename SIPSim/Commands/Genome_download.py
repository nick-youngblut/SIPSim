#!/usr/bin/env python

"""
genome_download: downloading genomes

Usage:
  genome_download [options] <accession_table>
  genome_download -h | --help
  genome_download --version

Options:
  <accessin_table>  Table of 'genome_name<tab>accession' entries.
                    No header. ('-' if from STDIN)
  -d=<d>            Output directory. [Default: .]
  -e=<e>            Email to use for NCBI queries. [Default: blank@gmail.com]
  -n=<n>            Number of cpus. [Default: 1]
  --debug           Debug mode (no multiprocessing).
  -h --help         Show this screen.
  --version         Show version.

Description:

"""

# import
## batteries
from docopt import docopt
from SIPSim import Genome_Download


# opt parse
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    Genome_Download.main(args)
   

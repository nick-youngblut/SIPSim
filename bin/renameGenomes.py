#!/usr/bin/env python

#--- Option parsing ---#
"""
renameGenomes: formatting genome sequences in a multi-fasta file for SIPsim

Usage:
  renameGenomes [options] <genome_fasta>
  renameGenomes -h | --help
  renameGenomes --version

Options:
  <genome_fasta>  Fasta file containing genome sequences ('-' if from STDIN).
  --quiet         Limit stderr output. 
  -h --help       Show this screen.
  --version       Show version.

Description:
  Reformating the sequence names so that they work with SIPSim.
  Edited fasta written to STDOUT.
  WARNING: there are no checks that the file is actually in fasta format.

"""

# import
## batteries
from docopt import docopt
import sys,os
import re


def main(Uargs):
    # input
    if Uargs['<genome_fasta>'] == '-':
        inFH = sys.stdin
    else:
        inFH = open(Uargs['<genome_fasta>'], 'rb')         

    # regexes
    re1 = re.compile(r'\W')
    re2 = re.compile(r'^_*(.*?)_*$')
    re3 = re.compile(r'_*complete_genome')
    re4 = re.compile(r'(.{78}).+')

    # iterating through file
    seq_names = dict()
    for line in inFH:
        line = line.rstrip()
        if line.startswith('>'):

            line = re1.sub('_', line)
            line = re2.sub(r'>\1', line)
            line = re3.sub('', line)
            line = re4.sub(r'\1', line)
            try:
                _ = seq_names[line]
                line = '_'.join(line, str(len(seq_names.keys())))
            except KeyError:
                pass

            seq_names[line] = 1
                
        print line

    # end
    inFH.close()
    
    
# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    main(Uargs)
    

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
import platform
import subprocess

## 3rd party
import parmap
import functools
import itertools

## application libraries
scriptDir = os.path.abspath(os.path.dirname(__file__))
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from Genome import Genome
import Utils


# functions
def get_os():
    """getting operating system; only unix-like machines"""
    
    OS = platform.uname()[0]
    if OS == 'Linux':
        OS = 'linux'
    elif OS == 'Darwin':
        OS = 'mac'
    else:
        sys.stderr.write('OS: "{}" not supported\n'.format(OS))

    return OS


def is_file(fileName):
    """Does file exist?"""
    if os.path.isfile(fileName) is False:
        raise IOError('"{}" does not exist'.format(fileName))

        
def sys_call(cmd, quiet=False):
    """System call of command."""
    try:
        if quiet:
            DEVNULL = open(os.devnull, 'w')
            proc = subprocess.Popen([cmd], shell=True, stdout=DEVNULL)
        else:
            proc = subprocess.Popen([cmd], shell=True)
    except subprocess.CalledProcessError:
        pass # handle errors in the called executable
    except OSError:
        raise OSError('No executable for command: "{}"\n'.format(cmd))

    output, err = proc.communicate()

    
def index_genome(taxonName, genomeFile, chilliDir, faToTwoBitExe, K_value, quiet=False):
    """indexing genome with MFEprimer indexing scripts. Just making system calls.
    Args:
    taxonName -- taxon name of genome fasta
    genomeFile -- file name of genome fasta
    chilliDir -- string with 'chilli' directory path
    faToTwoBitExe -- string with path of faToTwoBit
    K_value -- k value used for indexing
    quiet -- quiet all messages
    """
    # status
    if not quiet:
        sys.stderr.write('Indexing: "{}"\n'.format(taxonName))
        
    # begin indexing
    formatExe = os.path.join(chilliDir, 'UniFastaFormat.py')
    is_file(formatExe)
    cmd = '{} -i {}'.format(formatExe, genomeFile)
    sys_call(cmd)    
    
    # faToTwoBit
    is_file(faToTwoBitExe)
    cmd = '{} {} {}'.format(faToTwoBitExe,genomeFile + '.unifasta',
                            genomeFile + '.2bit')
    sys_call(cmd)
    
    # index
    indexExe = os.path.join(chilliDir, 'mfe_index_db.py')
    is_file(indexExe)
    cmd = '{} -f {} -k {}'.format(indexExe,
                                  genomeFile + '.unifasta',
                                  K_value)
    sys_call(cmd)

    # cleanup
    os.remove(genomeFile + '.unifasta')



def main(Uargs):
    # machine info
    OS = get_os()
    machine = platform.machine()
    
    # dirs
    chilliDir = os.path.join(libDir, 'chilli')
    if machine.endswith('_64'):
        faToTwoBitExe = os.path.join(libDir, 'bin', OS, '64/faToTwoBit')
    else:
        faToTwoBitExe = os.path.join(libDir, 'bin', OS, '32/faToTwoBit')        

        
    # loading genome list
    genomeList = Utils.parseGenomeList(Uargs['<genomeList>'], filePath=Uargs['--fp'])

    # setting function for parallel calling
    index_genome_par = functools.partial(index_genome,
                                         chilliDir=chilliDir,
                                         faToTwoBitExe=faToTwoBitExe,
                                         K_value=Uargs['--K_value'],
                                         quiet=Uargs['--quiet']
                                         )                                         

    # indexing genomes in parallel
    parmap.starmap(index_genome_par, genomeList, processes=int(Uargs['--np']), chunksize=1)

    # status
    print "#-- All genomes indexed --#"

    
    
# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    main(Uargs)
    

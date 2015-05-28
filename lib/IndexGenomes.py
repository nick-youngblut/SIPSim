#!/usr/bin/env python

# import
## batteries
import sys,os
import platform
## 3rd party
import functools
import dill
from pathos.multiprocessing import ProcessingPool
## application libraries
scriptDir = os.path.abspath(os.path.dirname(__file__))
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from Genome import Genome
import Utils

    
def index_genome(x, chilliDir, faToTwoBitExe, K_value, quiet=False):
    """indexing genome with MFEprimer indexing scripts. Just making system calls.
    Args:
    x -- (taxonName, genomeFile)
      taxonName -- taxon name of genome fasta
      genomeFile -- file name of genome fasta
    chilliDir -- string with 'chilli' directory path
    faToTwoBitExe -- string with path of faToTwoBit
    K_value -- k value used for indexing
    quiet -- quiet all messages
    """
    taxonName,genomeFile = x
    # status
    if not quiet:
        sys.stderr.write('Indexing: "{}"\n'.format(taxonName))
        
    # begin indexing
    formatExe = os.path.join(chilliDir, 'UniFastaFormat.py')
    Utils.is_file(formatExe)
    cmd = '{} -i {}'.format(formatExe, genomeFile)
    Utils.sys_call(cmd)    
    
    # faToTwoBit
    Utils.is_file(faToTwoBitExe)
    cmd = '{} {} {}'.format(faToTwoBitExe,genomeFile + '.unifasta',
                            genomeFile + '.2bit')
    Utils.sys_call(cmd)
    
    # index
    indexExe = os.path.join(chilliDir, 'mfe_index_db.py')
    Utils.is_file(indexExe)
    cmd = '{} -f {} -k {}'.format(indexExe,
                                  genomeFile + '.unifasta',
                                  K_value)
    Utils.sys_call(cmd)

    # cleanup
    os.remove(genomeFile + '.unifasta')



def main(Uargs):
    # machine info
    OS = Utils.get_os()
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
    pfunc = functools.partial(index_genome,
                              chilliDir=chilliDir,
                              faToTwoBitExe=faToTwoBitExe,
                              K_value=Uargs['--K_value'],
                              quiet=Uargs['--quiet'])                                         


    # indexing genomes in parallel
    pool = ProcessingPool(nodes=int(Uargs['--np']))
    if Uargs['--debug']:
        KDE_BD = map(pfunc, genomeList)
    else:
        KDE_BD = pool.map(pfunc, genomeList)

    # status
    sys.stderr.write('#-- All genomes indexed --#\n')

    
    
# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    main(Uargs)
    

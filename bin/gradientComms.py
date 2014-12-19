#!/usr/bin/env python

#--- Option parsing ---#
"""
gradientComms: simulate communities in the samples used for SIP

Usage:
  gradientComms [options] <genomeList>
  gradientComms -h | --help
  gradientComms --version

Options:
  <genomeList>  A file listing: taxonName<tab>genomeSeqFileName
  --fp=<fp>     Full path to genomeSeqFiles (if not in genomeList file).
  --pf=<pf>     Profile file setting optionsal arguments for grinder. See grinder help for more info.
  --exe=<exe>   Grinder executable. [Default: ./bin/grinder.pl]
  -h --help     Show this screen.
  --version     Show version.
  --debug       Debug mode

Description:
  Simple wrapper for calling modified version of grinder ('grinderSIP').
  grinderSIP just outputs a table of taxon abundances and does not actually
  simulate reads, so options that affect read simulation will have not affect.
  Options can be provided to grinderSIP with the profile file ('--pf').

  Output:
    A tab-delimited table of taxon abundances for each library is written to STDOUT.

"""

# import
## batteries
from docopt import docopt
import sys,os
import subprocess

## 3rd party
#import parmap

## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)



# functions
def call_grinder(genomeListFile, filePath=None, profileFile=None, exe='grinder'):
    """System call of grinder
    Args:
    genomeListFile -- genome list file name
    filePath -- full path to genome files (if needed)
    exe -- name of grinder executable
    profileFile -- profile file to provide to grinder
    """
    
    cmd = '{exe} -rf {glf} -fp {fp}'.format(exe=exe, glf=genomeListFile, fp=filePath)
    
    if profileFile is not None:
        cmd += ' -pf {}'.format(profileFile)

        
    #print os.system(cmd)
    print subprocess.check_output([cmd], shell=True)
    


# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')

    # grinder exe    
    grinderExe = str(args['--exe'])
    if grinderExe.startswith('./'):
        grinderExe = os.path.join(scriptDir, '.' + grinderExe)

        
    # call grinder; print output
    call_grinder(args['<genomeList>'],
                 filePath=args['--fp'],
                 profileFile=args['--pf'],
                 exe=grinderExe)
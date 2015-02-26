#!/usr/bin/env python

#--- Option parsing ---#
"""
fragGC: simulate genomeic fragments that would be found in a isopycnic gradient and calc frag G+C

Usage:
  fragGC [options] <genomeList>
  fragGC -h | --help
  fragGC --version

Options:
  <genomeList>  A file listing: taxonName<tab>genomeSeqFileName
  --fp=<fp>     Full path to genomeSeqFiles (if not in genomeList file).
  --rtr=<rtr>   Read template length range (min,max).
                [Default: 400,1200]
  --rtl=<rtl>   Read template length distribution (see Description).
                [Default: uniform,250,250]
  --nf=<nf>     Number of fragments to simulate per genome.
                If value ends in 'X', then the number of reads needed
                to obtain that coverage in the largest genome will be used
                (determined by simulating fragments using the specified
                fragment length distribution).
                [Default: 10000]
  --fld=<ld>    Fragment length distribution (see Description).
                [Default: skewed-normal,9000,2500,-5]
  --flr=<m>     Fragment length range (min,max). 
                [Default: 4000,None]
  --fr=<fr>     Fasta of forward & reverse primers (if amplicons).
  --np=<np>     Number of genomes to process in parallel.
                [Default: 1]  
  -h --help     Show this screen.
  --version     Show version.
  --debug       Debug mode

Description:
  Simulate the genomic fragments that would be found in an isopycnic gradient.
  The location and G+C of each simulated fragment is written to a table.

  Genome sequence file names in the <genomeList> file should either include 
  the full path to the file or the path needs to be provided by the '--fp' flag.

  '--rtr' and '--rtl':
    'read template' refers to the template nucleotide molecule from
    which the read originated. Template size is needed to constrain the location and lengths
    of the simulated genome fragments that contain these templates (G+C calculated from these fragments).
    Templates could either be an amplicon (e.g. 16S rRNA sequencing)
    or a genomic fragment (e.g. shotgun metagenomics). '--rtl' constrains
    the size of this template. If primers are provided, the in-silico generated amplicons
    are filtered to those that fit in the specified read template range.
    '--rtl' is used to determine the template lengths of shotgun metagenomic reads
    (not contrained by provided primers).

    ** Example for 16S rRNA amplicons **
    Assumning amplicon lengths of 500-650 bp, use: '--rtr 500,650'.
    This will filter out all in-silico amplicons outside of this range.


  **Distributions** 

  normal:
    Parameters: location (mean), scale (standard deviation)

  uniform:
    Parameters: low, high
    
  skewed-normal:
    Parameters: location (mean), scale (stardard deviation), shape (skew)
    Example: --fld  skewed-normal,11000,1000,-1000

  truncated-normal:
    Parameters: location, scale, low, high

  **Output**
    A tab-delim file written to STDOUT.

"""

# import
## batteries
from docopt import docopt
import sys,os
import re
import functools
import itertools

## 3rd party
import parmap

## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from Genome import Genome
from SimFrags import SimFrags
import Utils



# functions
def by_genome(inFile, taxonName, args):
    """All processing conducted per genome.
    Args:
    inFile -- genome sequence file name
    taxonName -- taxon name of genome
    args -- user-provided args as dict
    Return:
    2d-list -- for each fragment: [taxonName,scaf,start,end,GC]
    """
    # status
    sys.stderr.write('Processing: "{}"\n'.format(taxonName))
    
    # input check
    assert 'scriptDir' in args, '"scriptDir" not in args'
    
    
    # checking for MFEprimer.py executable
    MFEprimerExe = os.path.join(args['scriptDir'], 'MFEprimer.py')
    if not os.path.isfile(MFEprimerExe):
        raise IOError('Cannot find executable "{}"'.format(MFEprimerExe))


    # making genome object
    assert '--fr' in args, '"--fr" must be provided in args'
    genome = Genome(inFile, taxonName, args['--fr'])
    
    
    # sequenced read template location: amplicons
    if genome.get_primerFile() is not None:
        # in-silico PCR
        assert '--rtr' in args, '"--rtr" must be in args'
        genome.callMFEprimer(rtr=args['--rtr'], MFEprimerExe=MFEprimerExe)
    
        # filtering overlapping in-silico amplicons
        genome.filterOverlaps()
        
        # skip genome if no amplicons 
        if genome.get_nAmplicons == 0:
            return None

            
    # simulating fragments
    simFO = SimFrags(fld=args['--fld'], flr=args['--flr'], rtl=args['--rtl'])
    fragList = []
    for i in xrange(int(args['--nf'])):
        (scaf,fragStart,fragEnd,fragSeq) = simFO.simFrag(genome)

        fragList.append([genome.get_taxonName(),
                         scaf,
                         str(fragStart),
                         str(fragEnd),
                         str(genome.calcGC(fragSeq))
                     ])

    return fragList


def get_NFragsByCoverage(genomeList, coverage, fld, flr, rtl):
    """Determining the number of fragments to simulate
    based on the coverage of the largest genome in the provided
    list.
    Args:
    genomeList -- iterable of genome fasta file names
    coverage -- user-defined coverage
    fld -- fragment length distribution
    flr -- fragment length range
    rtl -- read template length
    Return:
    integer of number of fragments to simulate
    """
    # input
    coverage = float(coverage.rstrip('xX'))
    simFO = SimFrags(fld=fld, flr=flr, rtl=rtl)
    
    # making dict of genome lengths; finding max length
    gl = {taxonName : [Genome(inFile, taxonName).length, inFile]
          for inFile,taxonName in genomeList}
    gl_max_taxon = sorted(gl, key=gl.get, reverse=True)[0]
    gl_max_len = gl[gl_max_taxon][0]
    gl_max_file = gl[gl_max_taxon][1]

    # determining number of fragments needed based on simulated fragments
    fragLenCov = gl_max_len * coverage    
    genome = Genome(gl_max_file, gl_max_taxon)    
    fragLenTotal = 0
    nFrags = 0
    while 1:
        nFrags += 1
        (scaf,fragStart,fragEnd,fragSeq) = simFO.simFrag(genome)
        fragLenTotal += abs(fragEnd - fragStart + 1)
        if fragLenTotal >= fragLenCov:
            # status
            sys.stderr.write('Using user-defined genome coverage of the largest ' \
                             'genome to determine number of fragments to simulate.\n')
            sys.stderr.write('  Largest genome: {}\n'.format(gl_max_taxon))
            sys.stderr.write('  Genome size (bp): {}\n'.format(gl_max_len))    
            sys.stderr.write('  Number of fragments needed: {}\n'.format(nFrags))    

            # return
            return nFrags
    
        
# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')
    
    # splitting comma-separated args
    sp = re.compile(' *, *')    
    args = {k:sp.split(str(v)) if sp.search(str(v)) else v for k,v in args.items()}
    
    # adding to args
    args['scriptDir'] = libDir

    # list of genome files
    genomeList =  Utils.parseGenomeList(args['<genomeList>'], filePath=args['--fp'])

    # n-fragments: coverage or set number
    if args['--nf'].endswith('X') or args['--nf'].endswith('x'):
        args['--nf'] = get_NFragsByCoverage(genomeList, args['--nf'],
                                            fld=args['--fld'],
                                            flr=args['--flr'],
                                            rtl=args['--rtl'])
        
    # analyzing each genome (in parallel)    
    by_genome_part = functools.partial(by_genome, args=args)

    print '\t'.join(['taxon_name','scaffoldID','fragStart','fragEnd','GC'])            
    if args['--debug']:
        for x in itertools.starmap(by_genome_part,genomeList):
            pass
    else:
        ret = parmap.starmap(by_genome_part,
                             genomeList,
                             chunksize=1,
                             processes=int(args['--np']))

        for x in ret:
            for y in x:
                print '\t'.join(y)                
        
    

        
"""Fragment KDE classes"""

# import
## batteries
import sys,os
import math
import logging
import functools
import cPickle as pickle
## 3rd party
import numpy as np
import pandas as pd
import scipy.stats as stats
import mixture
import dill
from pathos.multiprocessing import ProcessingPool
# amplication
import SIPSimCython
from Genome import Genome
from SimFrags import SimFrags
import Utils

# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def load_frags_table(inFH, sep='\t'):
    """Loading frag info table as a dict of dicts of 2d lists.
    {taxon_name : {scaffold : [fragStart, fragEnd, GC]}}
    Args:
    inFH -- file handle
    sep -- value delimiter
    """    
    header_vals = set(['taxon_name','scaffoldID','fragStart',
                       'fragLength','fragGC'])
    
    d = dict()
    lineNum = 0
    for line in inFH.readlines():
        lineNum += 1
        line = line.rstrip().split(sep)

        #header
        if lineNum == 1:            
            if not (header_vals == set(line) or header_vals < set(line)):
                msg = 'The fragGC table does not have all'\
                      ' required columns:\n\t{}'\
                      .format(','.join(header_vals))
                raise IOError(msg)
            header_idx = {line[i]:i for i in xrange(len(line))}
        # body            
        else:
            taxon_name = line[header_idx['taxon_name']]
            try:
                type(d[taxon_name])
            except KeyError:
                d[taxon_name] = dict()
                d[taxon_name]['fragLength'] = []
                d[taxon_name]['fragGC'] = []

            fragLength = line[header_idx['fragLength']]
            fragGC = line[header_idx['fragGC']]
            d[taxon_name]['fragLength'].append(fragLength)
            d[taxon_name]['fragGC'].append(fragGC)
    return d

            
def load_frags_pickle(inFH):
    """Loading frag GC info assuming a pickled python object
    produced by SIPSim fragGC.
    Args:
    inFH -- file handle
    """
    fojb =  pickle.load(inFH)

    d = dict()
    for x in fojb:
        taxon_name = x[0]
        d[taxon_name] = dict()
        d[taxon_name]['fragLength'] = []
        d[taxon_name]['fragGC'] = []
            
        for scaf,v in x[1].items():            
            for z in v:
                # fragStart, fragLength, fragGC
                d[taxon_name]['fragLength'].append(z[1])
                d[taxon_name]['fragGC'].append(z[2])                
    return d



def load_frags(fileName):
    """Loading fragment data (pickled) table.
    Args:
    fileName -- name of the fragment data table
    Return:
    dict{dict} -- {taxon_name:{key:value}}
    """
    try:
        inFH = open(fileName, 'r')
    except IOError:
        inFH = sys.stdin

    try:
        frag_data = load_frags_pickle(inFH)
    except pickle.UnpicklingError:
        inFH.seek(0)
        frag_data = load_frags_table(inFH)            

    inFH.close()

    return frag_data


def fit_kde(frag_data, bw_method=None):
    """Returns multivariate KernelDensity function fit to
    fragment buoyant density (calculated from G+C) 
    and fragment lengths.
    Bandwidth selection based on bandwidth attribute.
    Args:
    frag_data -- dict of lists (fragment info)
    bw_method -- passed to stats.gaussian_kde
    Return:
    dict of kde objects {taxon_name:kde}
    """
    try:
        bw_method = float(bw_method)
    except TypeError:
        pass
    
    kdes = dict()
    for taxon_name,data in frag_data.items():
        # getting GC & length values
        try:
            frag_GC = data['fragGC']
        except KeyError:
            msg = 'Taxon: {}: cannot find "fragGC"'            
            raise KeyError, msg.format(taxon_name)
        try:
            frag_len = data['fragLength']
        except KeyError:
            msg = 'Taxon: {}: cannot find "fragLength"'            
            raise KeyError, msg.format(taxon_name)

        # GC2BD
        frag_BD = SIPSimCython.GC2BD(np.array(frag_GC))

        # kde fitting
        try:
            kdes[taxon_name] = stats.gaussian_kde([frag_BD, frag_len], 
                                                  bw_method=bw_method)
        except ValueError:
            kdes[taxon_name] = None

    return kdes


     
# functions
def by_genome(x, args):
    """All processing conducted per genome.
    Args:
    x -- [inFile,taxonName]
      inFile -- genome sequence file name
      taxonName -- taxon name of genome
    args -- user-provided args as dict
    Return:
    2d-list -- for each fragment: [taxonName,scaf,start,end,GC]
    """
    taxonName,inFile = x
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
    if genome.primerFile is not None:
        # in-silico PCR
        assert '--rtr' in args, '"--rtr" must be in args'
        genome.callMFEprimer(rtr=args['--rtr'], MFEprimerExe=MFEprimerExe)
    
        # filtering overlapping in-silico amplicons
        genome.filterOverlaps()
        
        
    # simulating fragments    
    simFO = SimFrags(fld=args['--fld'], flr=args['--flr'], rtl=args['--rtl'])
    nFragsMade = 0
    fragList = dict()
    ## if no amplicons
    if genome.nAmplicons == 0:
        pass
    ## if using coverage
    elif args['--nf'].endswith('X') or args['--nf'].endswith('x'):
        coverage = float(args['--nf'].rstrip('xX'))
        fragLenCov = genome.length * coverage
        fragLenTotal = 0
        while 1:
            (scaf,fragStart,fragLen,fragGC) = simFO.simFrag(genome)
            try:
                type(fragList[scaf])
            except KeyError:
                fragList[scaf] = []
                                
            if fragStart == "NA":
                break
            elif fragLenTotal > fragLenCov:
                break
            fragLenTotal += fragLen 

            nFragsMade += 1
            fragList[scaf].append([fragStart, fragLen, fragGC])            
    ## if using fixed number of fragments
    else:            
        for i in xrange(int(args['--nf'])):
            (scaf,fragStart,fragLen,fragGC) = simFO.simFrag(genome)

            try:
                type(fragList[scaf])
            except KeyError:
                fragList[scaf] = []

            if fragStart == "NA":
                break

            nFragsMade += 1
            fragList[scaf].append([fragStart, fragLen, fragGC])
                
    # status
    sys.stderr.write('  Genome name: {}\n'.format(genome.taxonName))                
    sys.stderr.write('  Genome length (bp): {}\n'.format(genome.length))
    if args['--nf']:
        msg = '  Number of amplicons: {}\n'
        sys.stderr.write(msg.format(genome.nAmplicons))
    msg = '  Number of fragments simulated: {}\n'
    sys.stderr.write(msg.format(nFragsMade))
                
    return [genome.taxonName, fragList]


def write_fragList(fragList):
    """Writing out fragList as a tab-delim table.
    """
    print '\t'.join(['taxon_name','scaffoldID','fragStart',
                     'fragLength','fragGC'])            
    for x in fragList:
        taxon_name = x[0]
        for scaf,v in x[1].items():
            for y in v:                
                print '\t'.join([taxon_name, scaf] + [str(i) for i in y])


def main(args):
    # list of genome files
    genomeList =  Utils.parseGenomeList(args['<genomeList>'], 
                                        filePath=args['--fp'])
        
    # analyzing each genome (in parallel)    
    pfunc = functools.partial(by_genome, args=args)
    
    # difussion calc in parallel
    pool = ProcessingPool(nodes=int(args['--np']))
    if args['--debug']:
        fragList = map(pfunc, genomeList)
    else:
        fragList = pool.map(pfunc, genomeList)

    # writing out table
    if args['--tbl']:
        write_fragList(fragList)
    else:
        dill.dump(fragList, sys.stdout)
        


"""Fragment KDE classes"""


#########
# Depreciated
#######

# import
## batteries
import sys
import math
import logging
import cPickle as pickle
## 3rd party
import numpy as np
import pandas as pd
import scipy.stats as stats
import mixture

# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def load_fragGC_table(inFH, sep='\t'):
    """Loading fragGC table as a dict of dicts of 2d lists.

    Parameters
    ----------
    inFH : file handle
    sep : str
        value delimiter

    Returns
    -------
    fragGC -- fragGC object
        {taxon_name : {scaffold : [fragStart, fragEnd, GC]}}
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

            
def load_fragGC_pickle(inFH):
    """Loading fragGC info assuming a pickled python object
    produced by SIPSim fragGC.
    
    Parameters
    ----------
    inFH : file handle
    
    Returns
    -------
    fragGC object
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


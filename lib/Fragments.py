"""Fragment KDE classes"""

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
#from sklearn.neighbors.kde import KernelDensity
import mixture
# amplication
import SIPSimCython

# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def load_frags_table(inFH, sep='\t'):
    """Loading frag info table as a dict of dicts of 2d lists.
    {taxon_name : {scaffold : [fragStart, fragEnd, GC]}}
    Args:
    inFH -- file handle
    sep -- value delimiter
    """    
    header_vals = set(['taxon_name','scaffoldID','fragStart','fragLength','fragGC'])
    
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


     
class Frag_multiKDE(object):
    """Multivariate KDE fit to fragment G+C and fragment lengths.
    Method:
    * load all fragGC data for >= taxon
    * fit mulit-var KDE to distribution (using gaussian_kde from scipy.stats)
    """
    
    def __init__(self, fragGC_file, bandwidth=None):
        """Using sklearn.neighbors.kde for distribution fitting. See
        the doc for that class to get more info on the parameters passed
        to it. 
        Args:
        fragGC_file -- File produced by fragGC subcommand
        bandwidth -- KDE bandwidth parameter passed to KDE function
        """
        self.fragGC_file = fragGC_file
        self.bandwidth = bandwidth

        # loading fragment values        
        def _load_fragGC(inFH):
            try:
                self._frag_data = load_fragGC_pickle(inFH)
            except pickle.UnpicklingError:
                inFH.seek(0)
                self._frag_data = load_fragGC_table(inFH)            
        try:
            with open(self.fragGC_file, 'r') as inFH:
                _load_fragGC(inFH)
        except IOError:
            _load_fragGC(sys.stdin)

        # iter by taxon to fit kde
        self._fragKDE = dict()
        for taxon_name in self._frag_data.keys():
            self._fragKDE[taxon_name] = self._fitKDE(taxon_name)

        
    def __repr__(self):
        out = ''        
        for taxon_name in self.iter_taxon_names():
            for kde in self.iter_kde([taxon_name]):
                out += '{}\n\t{}\n'.format(taxon_name, kde.__repr__())
        return out


    def iter_taxon_names(self):
        """Iterate through all taxon names.
        """
        for taxon_name in self._fragKDE.keys():
            yield taxon_name

            
    def iter_kde(self, taxon_names):
        """Getting the kde functions for all taxa-listed.
        Args:
        taxon_names -- list of taxon names
        """
        assert not isinstance(taxon_names, str), 'Provide a list-like object'

        for taxon_name in taxon_names:
            try:
                yield self._fragKDE[taxon_name]
            except KeyError:
                raise KeyError('No KDE fit for taxon "{}"'.format(taxon_name))

    
    def get_all_kde(self):
        """Return all KDE objects as list [[taxon_name, kde],]
        """
        return [[taxon_name, self._fragKDE[taxon_name]] 
                 for taxon_name in self.iter_taxon_names()]

                
    def sampleTaxonKDE(self, taxon_name, size=1):
        """Sample from the KDE function set for a taxon.
        Args:
        taxon_name -- name of taxon
        size -- number of fragments to sample
        #args and kwargs -- passed to KDE function
        Return:
        #2d numpy array -- [[frag_GC,frag_len], ...]
        generator: 1d numpy array [frag_GC, frag_len]
        """
        # asserting kde function
        try:            
            kde = self._fragKDE[taxon_name]
        except KeyError:
            raise KeyError('taxon "{}" does not have a KDE function set'.format(taxon_name))

        if kde is None:
            return None
        else:
            return kde.resample(size=size)
                
                                   
    def _fitKDE(self, taxon_name):
        """Returns multivariate KernelDensity function fit to fragment GC
        values and lengths.
        Bandwidth selection based on bandwidth attribute.
        Args:
        taxon_name -- taxon name string
        """
        #vals = np.vstack([taxon_df.GC, abs(taxon_df.fragEnd - taxon_df.fragStart)])        
        vals = [self._get_fragGC(taxon_name),
                self._get_fragLength(taxon_name)]

        try:
            return stats.gaussian_kde(vals, bw_method=self.bandwidth)
        except ValueError:
            return None

        
    def _iter_df_by_taxon_name(self):
        """Iterate by unique taxon names.
        """
        taxon_names = self.fragGC_df['taxon_name'].unique()
        for taxon_name in taxon_names:
            yield (taxon_name, self.fragGC_df.loc[self.fragGC_df['taxon_name'] == taxon_name])


    def _get_fragGC(self, taxon_name):
        """Getting fragGC values for a taxon.
        """
        assert hasattr(self, '_frag_data'), "No _frag_data attribute"
        try:
            type(self._frag_data[taxon_name])
        except KeyError:
            raise KeyError('taxon: "{}" not found'.format(taxon_name))
        
        return self._frag_data[taxon_name]['fragGC']
            

    def _get_fragLength(self, taxon_name):
        """Getting fragGC values for a taxon.
        """
        assert hasattr(self, '_frag_data'), "No _frag_data attribute"
        try:
            type(self._frag_data[taxon_name])
        except KeyError:
            raise KeyError('taxon: "{}" not found'.format(taxon_name))
        
        return self._frag_data[taxon_name]['fragLength']


        

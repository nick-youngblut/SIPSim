"""Fragment KDE classes"""

# import
## batteries
import sys
import math
import logging 
## 3rd party
import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn.neighbors.kde import KernelDensity
import pymix.mixture as mixture

# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


class Frag_multiKDE(object):
    """Multivariate KDE fit to fragment G+C and fragment lengths.
    Method:
    * load all fragGC values for a taxon
    * fit mulit-var KDE to distribution (using gaussian_kde from scipy.stats)
    """
    
    def __init__(self, fragGC_file, bandwidth=0.0001):
        """Using sklearn.neighbors.kde for distribution fitting. See
        the doc for that class to get more info on the parameters passed
        to it. 
        Args:
        fragGC_file -- File produced by fragGC subcommand
        bandwidth -- KDE bandwidth parameter passed to KDE function
        """
        self.fragGC_file = fragGC_file
        self.bandwidth = bandwidth

        # loading values as pandas dataframe
        self.fragGC_df = pd.read_csv(self.fragGC_file, sep='\t')
        
        # iter by taxon to fit kde
        self._fragKDE = dict()
        for (taxon_name,taxon_df) in self._iter_df_by_taxon_name():
            self._fragKDE[taxon_name] = self._fitKDE(taxon_df)

        # delete df
        self.fragGC_df = None

        
    def __repr__(self):
        out = ''        
        for taxon_name in self.iter_taxon_names():
            for kde in self.iter_fragGC_kde([taxon_name]):
                out += '{}\n\t{}\n'.format(taxon_name, kde.__repr__)
        return out


    def iter_taxon_names(self):
        """Iterate through all taxon names.
        """
        for taxon_name in self._fragGC_kde.keys():
            yield taxon_name

            
    def iter_fragGC_kde(self, taxon_names):
        """Getting the kde functions for all taxa-listed.
        Args:
        taxon_names -- list of taxon names
        """
        assert not isinstance(taxon_names, str), 'Provide a list-like object'

        for taxon_name in taxon_names:
            try:
                yield self._fragGC_kde[taxon_name]
            except KeyError:
                raise KeyError('No KDE fit for taxon "{}"'.format(taxon_name))

                
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

#        print kde.resample(size=10).T; sys.exit()
            
        # return samples
        #return kde.resample(*args, **kwargs).T
        # yeild samples
        #return (kde.resample(size=1)[:,0] for x in xrange(size))
        return kde.resample(size=size)
                
                                   
    def _fitKDE(self, taxon_df):
        """Returns multivariate KernelDensity function fit to fragment GC values and lengths.
        #Also sets 'kde_element_len', which signifies the number of values returned
        #as 1 'sample' from the *.sample() method. This is needed because 1 'sample'
        #can actually return many (100's, 1000's or more) values depending on the
        #values used to fit the KDE.
        """
        vals = np.vstack([taxon_df.GC, abs(taxon_df.fragEnd - taxon_df.fragStart)])                                 
        return stats.gaussian_kde(vals, bw_method=self.bandwidth)

        
    def _iter_df_by_taxon_name(self):
        """Iterate by unique taxon names.
        """
        taxon_names = self.fragGC_df['taxon_name'].unique()
        for taxon_name in taxon_names:
            yield (taxon_name, self.fragGC_df.loc[self.fragGC_df['taxon_name'] == taxon_name])

        

class FragGC_KDE(object):
    """Distributions fit to the fragment G+C values produced by the fragGC subcommand.
    Method:
    * load all fragGC values for a taxon
    * fit KDE to distribution (using KDE from sklearn.neighbors)
    """
    
    def __init__(self, fragGC_file, kernel='gaussian', bandwidth=0.0001):
        """Using sklearn.neighbors.kde for distribution fitting. See
        the doc for that class to get more info on the parameters passed
        to it. 
        Args:
        fragGC_file -- File produced by fragGC subcommand
        bandwidth -- KDE bandwidth parameter passed to KDE function
        kernel -- kernal parameter passed to KDE function
        """
        self.fragGC_file = fragGC_file
        self.kernel = kernel
        self.bandwidth = bandwidth

        # loading values as pandas dataframe
        self.fragGC_df = pd.read_csv(self.fragGC_file, sep='\t')
        
        # iter by taxon to fit kde
        self._fragGC_kde = dict()
        for (taxon_name,taxon_df) in self._iter_df_by_taxon_name():
            self._fragGC_kde[taxon_name] = self._fitDist2GC(taxon_df)

        # delete df
        self.fragGC_df = None

                        
    def __repr__(self):
        out = ''        
        for taxon_name in self.iter_taxon_names():
            for kde in self.iter_fragGC_kde([taxon_name]):
                out += '{}\n\t{}\n'.format(taxon_name, kde.__repr__)
        return out
        
        
    def iter_taxon_names(self):
        """Iterate through all taxon names.
        """
        for taxon_name in self._fragGC_kde.keys():
            yield taxon_name

            
    def iter_fragGC_kde(self, taxon_names):
        """Getting the kde functions for all taxa-listed.
        Args:
        taxon_names -- list of taxon names
        """
        assert not isinstance(taxon_names, str), 'Provide a list-like object'

        for taxon_name in taxon_names:
            try:
                yield self._fragGC_kde[taxon_name]
            except KeyError:
                raise KeyError('No KDE fit for taxon "{}"'.format(taxon_name))

                
    def sampleTaxonKDE(self, taxon_name, *args, **kwargs):
        """Sample from the KDE function set for a taxon.
        Args:
        taxon_name -- name of taxon
        args and kwargs -- passed to KDE function
        Return:
        list -- sampled values
        """
        # edit n_samples by how many samples returned from 1 'n_sample'
        try:
            n_samples = kwargs['n_samples']
            del(kwargs['n_samples'])
        except KeyError:
            n_samples = 1

        n_samples_act = int(math.ceil(n_samples / self.kde_element_len))
        if n_samples_act < 1: n_samples_act = 1

        # calling kde.sample()
        try:            
            kde = self._fragGC_kde[taxon_name]
        except KeyError:
            raise KeyError('taxon "{}" does not have a KDE function set'.format(taxon_name))

        # return
        return kde.sample(*args, n_samples=n_samples_act, **kwargs).flatten()[0:n_samples]
                
                                   
    def _fitDist2GC(self, taxon_df):
        """Returns KernelDensity function fit to GC values.
        Also sets 'kde_element_len', which signifies the number of values returned
        as 1 'sample' from the *.sample() method. This is needed because 1 'sample'
        can actually return many (100's, 1000's or more) values depending on the
        values used to fit the KDE.
        """
        self.kde_element_len = len(taxon_df['GC'])
        return KernelDensity(kernel=self.kernel, bandwidth=self.bandwidth).fit(taxon_df['GC'])

        
    def _iter_df_by_taxon_name(self):
        """Iterate by unique taxon names.
        """
        taxon_names = self.fragGC_df['taxon_name'].unique()
        for taxon_name in taxon_names:
            yield (taxon_name, self.fragGC_df.loc[self.fragGC_df['taxon_name'] == taxon_name])

            

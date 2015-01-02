# import
## batteries
import os,sys
import functools
import itertools
import math
import re
import logging
from collections import defaultdict
## 3rd party
import pandas as pd
import numpy as np
import scipy.stats as stats
from sklearn.neighbors.kde import KernelDensity
import pymix.mixture as mixture
from intervaltree import Interval, IntervalTree


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

                
    def sampleTaxonKDE(self, taxon_name, *args, **kwargs):
        """Sample from the KDE function set for a taxon.
        Args:
        taxon_name -- name of taxon
        args and kwargs -- passed to KDE function
        Return:
        2d numpy array -- [[frag_GC,frag_len], ...] 
        """

        # asserting kde function
        try:            
            kde = self._fragKDE[taxon_name]
        except KeyError:
            raise KeyError('taxon "{}" does not have a KDE function set'.format(taxon_name))

        # return samples
        return kde.resample(*args, **kwargs).T
                
                                   
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

            

class _table(object):
    """Template class for reading in SIPSim tables"""
    
    def __init__(self, df, filename):
        """
        Args:
        df -- pandas dataframe
        """
        self.df = df
        self.tableFileName = filename

        # library as string
        try:
            self.df['library'] = self.df['library'].astype(str)
        except KeyError:
            raise KeyError('"library" column not found in table: "{}"!'.format(filename))

            
    # load from csv
    @classmethod
    def from_csv(cls, filename, **kwargs):
        """Read in table file to a pandas dataframe.
        Args:
        filename -- Table file name
        kwargs -- passed to pandas.read_csv
        """
        df = pd.read_csv(filename, **kwargs)
        return cls(df, filename)

        
    # get/set/iter
    def iter_uniqueColumnValues(self, columnID):
        """General iteration of unique column values.
        """
        try:
            for l in self.df[columnID].unique():
                yield l
        except KeyError:
            raise KeyError('Column "{}" not found'.format(columnID))
            
    def iter_libraries(self):
        """iterating through all unique library IDs
        """
        for libID in self.iter_uniqueColumnValues('library'):
            yield libID
                
    def iter_taxa(self, libID=None):
        """Iterating through all unique taxon names.
        """
        if libID is None:
            for taxon_name in self.iter_uniqueColumnValues('taxon_name'):
                yield taxon_name
        else:
            df_lib = self.df.loc[self.df['library'] == libID]
            for taxon_name in df_lib['taxon_name'].unique():
                yield taxon_name
                
    def iter_taxonRowsInLib(self, libID):
        """Iterate through all subset dataframes containing just 1 taxon.
        Args:
        libID -- str; library ID        
        """
        df_lib = self.df.loc[self.df['library'] == libID]
        for taxon_name in df_lib['taxon_name'].unique():
            yield (taxon_name, df_lib.loc[df_lib['taxon_name'] == taxon_name])
                
    def __repr__(self):
        return self.df.__repr__()


        
class CommTable(_table):
    
    def __init__(self, *args, **kwargs):
        _table.__init__(self, *args, **kwargs)

        
    def set_abs_abund(self, abs_abund):
        """Setting the absolute abundance of each taxon based on user-value.
        Args:
        abs_abund -- int; total abundance of all taxa in each library
        """
        try:
            self.df['abs_abund'] = self.df['rel_abund_perc'] / 100 * int(abs_abund)
        except KeyError:
            raise KeyError('"rel_abund_perc" column not found in comm file')


    def taxonInLib(self, taxon_name, libID):
        """Checking whether taxon in the selected library
        Args:
        taxon_name -- taxon name
        libID -- library ID
        Return:
        boolean
        """
        libID = re.sub('\D+', '', libID)
        dfSub = self.df.loc[self.df['taxon_name'] == taxon_name]
        libs = [str(x) for x in dfSub['library'].tolist()]
        return libID in libs

            
    def get_taxonAbund(self, libID, taxon_name, rel=False):
        """Getting the abundance of taxon in the comm file.
        Absolute abundance is returned unless rel=True.
        If rel=True, relative abundance is returned.
        Args:
        rel -- return relative abundance instead of absolute?
        """
        if rel == True:
            retCol = 'rel_abund_perc'
        else:
            retCol = 'abs_abund'
            
        assert retCol in self.df.columns, \
            '"{}" column not found'.format(retCol)
        df_sub =  self.df.loc[(self.df['library'] == libID) &
                              (self.df['taxon_name'] == taxon_name)]
        return int(df_sub[retCol])

    
    def get_unique_taxon_names(self):
        return self.df['taxon_name'].unique()        

        
class IsoIncorpTable(_table):
    """Isotope incorporation class; provided isotope table is converted to pymix functions"""

    def __init__(self, *args, **kwargs):
        _table.__init__(self, *args, **kwargs)

        # status
        self._weightWarn = []
        
        # setting column types
        self.df['weight'] = self.df['weight'].astype(float)
        self.df['param_value'] = self.df['param_value'].astype(float)

        # setting incorp functions
        self.incorpFuncs = self._initialize()
        
        
    def _initialize(self):
        """For each library-taxon, making distribution function instance.
        Return:
        dict -- libID : taxon_name : distribution_function
        """
        incorpFuncs = defaultdict(dict)
        for libID in self.iter_libraries():
            for (taxon_name,taxon_df) in self.iter_taxonRowsInLib(libID):
                incorpFuncs[libID][taxon_name] = self._get_incorpDistFunc(taxon_df, libID, taxon_name)

        return incorpFuncs
    
        
    def _get_incorpDistFunc(self, df, libID, taxon_name):
        """Based on dataframe of parameters for a taxon's isotope incorporation
        (intra-population incorporation), assign a pymix distribution function.
        Args:
        df -- pandas dataframe; just for 1 taxon in 'taxon_name' column
        libID -- str; library ID
        taxon_name -- str; taxon name
        Return:
        scipy mixture model function
        """
        assert len(df['taxon_name'].unique()) == 1, '>=1 taxon name provided'

        # supported standard distribution functions
        psblFuncs = {'normal' : mixture.NormalDistribution,
                     'uniform' : mixture.UniformDistribution}
        
        # storing pymix distribution functions & weights for each distribution
        allDistFuncs = []
        allDistWeights = []
        
        # iter by distribution
        allDists = df[['distribution','distribution_index']].drop_duplicates()        
        for (distID, distIdx) in zip(allDists['distribution'], allDists['distribution_index']):
            df_sub = df.loc[df['distribution_index'] == distIdx]

            params = {k:v for k,v in
                      zip(df_sub['param'], df_sub['param_value'])}

            # rounding may cause intra-pop start > end, but should be ~exact
            try: 
                startParam = params['start']
                endParam = params['end']
            except KeyError:
                pass
            else:
                if startParam > endParam:
                    params['start'] = endParam
                    params['end'] = startParam
                elif startParam == endParam:
                    if startParam >= 100:
                        params['start'] -= 1e-10
                    else:
                        params['end'] += 1e-10

            # getting distribution
            try:
                distFunc = psblFuncs[distID.lower()]
                allDistFuncs.append(distFunc(**params))
            except KeyError:
                msg = 'Distribution "{}" not supported'
                raise KeyError(msg.format(distID))
            except TypeError:
                msg = 'For: {}=>{}, distribution "{}" does not support all provided params: "{}"'
                raise TypeError(msg.format(libID, taxon_name, distID, ','.join(params.keys())))
                
            # getting weights
            weights = [x for x in df_sub['weight'].unique()]
            assert len(weights) == 1, \
                '>=1 weight set for: {}=>{}=>{}'.format(libID, taxon_name, distID)
            allDistWeights.append(weights[0])

        # all weights == 1?
        if sum(allDistWeights) != 1:
            if libID in self._weightWarn:
                pass
            else:
                msg = 'lib->{}: intra-population weights do not sum to 1.' +\
                      ' Setting equal weights.'
                logging.warning(msg.format(libID))
                self._weightWarn.append(libID)
            
            allDistWeights = [1.0 / len(allDistWeights)] * len(allDistWeights)
            
        assert len(allDistWeights) == len(allDistFuncs), \
            'len(weights) != len(distFuncs) for: {}=>{}'.format(libID, taxon_name)
        
        # returning mixture model
        return mixture.MixtureModel(len(allDistFuncs), allDistWeights, allDistFuncs)
        
        
    def iter_incorpFuncs(self):
        """Iterating through all incorporation distribution functions.
        Yields:
        libraryID, taxon_name, distribution_function
        """
        for (libID,v1) in self.incorpFuncs.items():
            for (taxon_name,v2) in v1.items():
                yield libID, taxon_name, v2

                
    def sample_incorpFunc(self, libID, taxon_name, n_samples=1):
        """Sampling from the intra-population isotope incorporation function
        set for the user-selected library-taxon.
        Args:
        libID -- library ID
        taxon_name -- name of taxon in the library
        n_samples -- number of samples to pull from distribution
        Return:
        iterator of incorporation percentages
        """
        try:
            incorpFunc = self.incorpFuncs[libID][taxon_name]
        except KeyError:
            raise KeyError('Cannot find library-taxon: "{}"-"{}"'.format(libID, taxon_name))

        for i in xrange(n_samples):
            yield incorpFunc.sample()[0]
                

class FracTable(_table):
    """class for fraction table (output of fractions subcommand)"""
    
    def __init__(self, *args, **kwargs):
        _table.__init__(self, *args, **kwargs)
        
        # 'delete' inherited methods involving taxa
        iter_taxa = property()
        iter_taxonRowsInLib = property()

        # assertions on columns
        assertCols = ['library','fraction','BD_min','BD_max','fraction_size']
        for col in assertCols:
            assert col in self.df.columns, \
                'Column "{}" not found in file "{}"'.format(col, self.tableFileName)

        # fraction ID as string
        self.df['fraction'] = self.df['fraction'].astype(str)
            
        # creating an interval tree of fractions
        self._set_itrees()

        
    def _set_itrees(self):
        """Setting a dict of interval trees for the BD-min/max ranges (1 per lib).
        Max in interval ranges is non-inclusive.
        """
        self.itrees = dict()

        def make_frac_interval(x):
            return Interval(float(x[0]),float(x[1]),x[2])
        
        for libID in self.iter_libraries():
            df_sub = self.df.loc[self.df['library'] == libID, ['BD_min','BD_max','fraction']]
            self.itrees[libID] = IntervalTree(list(df_sub.apply(make_frac_interval, axis=1)))

            
    def which_frac(self, libID, BD_value):
        """Determine which of the simulated fractions that the BD value falls into.
        Using an interval tree of BD-min/max values.
        BD_max is non-inclusive.
        Args:
        libID -- library ID
        BD_value -- Bouyant density value
        Return:
        list of fractionIDs that the BD_value falls into
        """
        libID = str(libID)
        BD_value = float(BD_value)
        
        assert hasattr(self, 'itrees'), 'Cannot find itrees attrib'
        
        try:
            fracIDs = [x.data for x in self.itrees[libID][BD_value]]            
        except KeyError:
            raise KeyError('Library "{}" not found in fraction table'.format(libID))

        fracIDs_len = len(fracIDs)
        if fracIDs_len == 0:
            return None
        elif fracIDs_len == 1:
            return fracIDs[0]
        else:
            msg = 'BD value "{}" in library "{}" falls into >1 fraction!'
            raise ValueError(msg.format(str(BD_value), libID))

                    
    def fracID2BDminmax(self, libID, fracIDs):
        """Convert fractionIDs to 'BDmin-BDmax
        Args:
        libID -- str; library ID
        fracIDs -- iterable; fraction IDs
        """
        df_sub = self.df.loc[(self.df['library'] == libID)]
        df_sub.index = df_sub['fraction']
        df_sub = df_sub.reindex(fracIDs)

        
        def catCol(x, sep='-', cols=['BD_min','BD_max']):
            return sep.join([libID] + [str(y) for y in x[cols]])
        
        return  list(df_sub.apply(catCol, axis=1))

        
    def get_fracBDminxax(self, libID, fracID):
        """Getting the BD min/max values for a fraction in a library.
        Args:
        libID -- str; library ID
        fracID -- str; fraction ID
        Return:
        list -- [BD_min, BD_max]
        """
        df_sub = self.df.loc[(self.df['library'] == libID) & (self.df['fraction'] == fracID)]
        return list(df_sub.loc[0][['BD_min','BD_max']])


    def get_fracIDs(self):
        """Return list of library-fractionIDs.
        Format: [libraryID]_[fractionID]_[BD_min]_[BD_max]
        """
        df_sub = self.df[['library','fraction','BD_min','BD_max']]
        return list(df_sub.apply(lambda x : '_'.join([str(y) for y in x]), axis=1))

    
    def get_libFracIDs(self, libID):
        """Return list of fractionIDs in a library
        Args:
        libID -- library ID
        """
        df_sub = self.df.loc[self.df['library'] == libID]
        return list(df_sub['fraction'])
        
        
    def get_libFracDict(self):
        """Return {libraryID : {fractionID : dict}}
        """
        libFrac = defaultdict(defaultdict(dict))
        for i in xrange(self.df.shape[0]):
            row = self.df.loc[i]
            libFrac[str(row['library'])][str(row['fraction'])] = dict()
        return libFrac
        
                
        
class OTU_table(object):
    """Class for the simulated OTU table"""
    
    def __init__(self, frac, g_noise, scale, abund_weight, isotope):
        """
        Args:
        frac -- Fractions class instance
        g_noise -- str; gradient noise distribution name (from scipy.stats)
        scale -- float; scale parameter for gradient noise distribution
        abund_weight -- float; isotope incorp abundance-weighting factor
        isotope -- str; name of isotope
        """
        self.g_noise = g_noise
        self.scale = scale
        self.abund_weight = abund_weight
        self.isotope = isotope.upper()
        
        # setting gradient 'noise' function
        self.g_noise_func = self._set_g_noise_func(self.g_noise, self.scale)
        
        # set isotope theoretical max BD
        self.isotopeMaxBD = self._set_isotopeMaxBD(self.isotope)

        
    def _set_g_noise_func(self, dist_name, scale):
        """Setting the gradient noise function as scipy distribution function.
        Args:
        dist_name -- name of scipy distribution function
        scale -- scale parameter of distribution
        """
        scale = float(scale)
        
        psblFuncs = {'cauchy' : stats.cauchy,
                     'normal' : stats.norm,
                     'uniform' : stats.uniform}

        try:
            func = psblFuncs[dist_name]
        except KeyError:
            raise KeyError('Distribution "{}" not supported.'.format(dist_name))

        return functools.partial(func, scale=scale)
                     
                     
    def _set_isotopeMaxBD(self, isotope):
        """Setting the theoretical max BD shift of an isotope (if 100% incorporation).
        Args:
        isotope -- str; name of isotope
        """
        psblIsotopes = {'13C' : 0.036,
                        '15N' : 0.016}
        try:
            return psblIsotopes[isotope]
        except KeyError:
            raise KeyError('Isotope "{}" not supported.'.format(isotope))
            

    def sample_g_noise_func(self, loc, n_samples=1):
        """Sampling the gradient noise distribution function.
        Args:
        loc -- loc param of the scipy.stats distribution function
        n_samples -- number of samples to draw from distrubtion
        Return:
        list -- values
        """
        return self.g_noise_func(loc=float(loc)).rvs(n_samples)
            
            
    def checkLibOverlap(self, libList):
        """Checking that libraries in dataframes fully overlap
        Args:
        libList -- list of list of unique libs (each lib list is compared to the others)
        """        
        return len(set(itertools.chain(*libList))) == len(set(libList[0])) 

        
    def make_emptyCountTable(self, taxon_names, fractionIDs):
        """Making a pandas dataframe (taxon x fraction) of zeros.
        Args:
        taxon_names -- iterable; all taxon names
        fractionIDs -- iterable; all fraction IDs
        """
        shape = (len(taxon_names),len(fractionIDs))
        return pd.DataFrame(np.zeros(shape), columns=fractionIDs, index=taxon_names)
        
        
    def get_isotopeMaxBD(self):
        return self.isotopeMaxBD
        
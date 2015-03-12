# import
## batteries
import os, sys
import pandas as pd
import numpy as np
from functools import partial
from collections import Counter
## application
from Utils import _table


class OTU_table(_table):

    def __init__(self, *args, **kwargs):
        _table.__init__(self, *args, **kwargs)


    def set_samp_dist(self, samp_dist, samp_dist_params):
        """Setting subsampling size distribution & params.
        Args:
        samp_dist -- str of numpy.random distribution attribute
        samp_dist_params -- dict of params for the distribution function
        """
        self.samp_dist_params = samp_dist_params        
        self.samp_dist = samp_dist


    def get_comm_size_stats(self):
        """Getting stats on the size of each community.
        Return: [min, mean, median, max]
        """
        counts = self.df.groupby(['library','fractions']).sum()['count']
        return [np.min(counts), np.mean(counts), np.median(counts), np.max(counts)]

        
    def subsample(self, no_replace=False):
        """Subsampling from each community.
        Using numpy.random.choice with taxon abundances as weights
        Args:
        no_replace -- subsample without replacement
        Return:
        pandas DataFrame of subsampled community
        """
        # assertions
        assert hasattr(self, 'samp_dist'), 'samp_dist attribute not found'
        assert hasattr(self, 'samp_dist_params'), 'samp_dist_params attribute not found'

        # all taxa
        all_taxa = Counter({x:0 for x in self.iter_taxa()})
        
        # subsampling
        df_sub = pd.DataFrame(columns=self.df.columns)
        for libID in self.iter_libraries():
            for fracID in self.iter_fractions(libID=libID):
                # single community
                comm = self.get_comm(libID, fracID)
                sub_comm = None

                # if no counts for community
                if np.sum(comm['count']) < 1:
                    sub_comm = comm.copy()
                else:
                    # size to subsample
                    samp_size = self.samp_dist(size=1)
                    
                    # sampling
                    counts = comm['count']
                    try:
                        sub_comm = Counter(np.random.choice(comm['taxon'],
                                                         size=samp_size,
                                                         replace= not no_replace,
                                                         p=counts/np.sum(counts)))

                        
                        # setting all taxa in counts
                        sub_comm.update(all_taxa)
                        
                        # count to dataframe
                        sub_comm = pd.DataFrame(sub_comm.items())
                        sub_comm.columns = ['taxon','count']
                        sub_comm.loc[:,'library'] = libID
                        sub_comm.loc[:,'fractions'] = fracID
                    except ValueError:
                        sub_comm = comm.copy()
                        comm.loc[:,'count'] = 0

                df_sub = pd.concat([df_sub, sub_comm])

                        

        df_sub['count'] = df_sub['count'].astype(int)
        return df_sub.reindex_axis(['library','fractions','taxon','count'], axis=1)
                
        
    def _same_low_high(self):
        """Returns low/high value if samp_dist_params contain keys 'low' and 'high
        and the values of these keys are equal.
        Else, returns False.
        """
        try:
            if self.samp_dist_params['high'] == self.samp_dist_params['low']:
                return self.samp_dist_params['high']
        except KeyError:
            return False
        except AttributeError:
            return False


    def get_comm(self, libID, fracID):        
        return self.df.loc[(self.df['library'] == libID) &
                           (self.df['fractions'] == fracID),:]
            
    # iter
    def iter_fractions(self, libID):
        """iterating through fractions of a certain library.
        Yields a fraction ID.
        """
        df_sub = self.df.loc[self.df['library'] == libID]
        for fracID in df_sub['fractions'].unique():
            yield fracID
        

    
    # properties/setters
    @property
    def samp_dist_params(self):
        return self._samp_dist_params
    @samp_dist_params.setter
    def samp_dist_params(self, x):
        self._samp_dist_params = x
        
    @property
    def samp_dist(self):
        return self._samp_dist
    @samp_dist.setter
    def samp_dist(self, x):
        self.samp_dist_name = x

        # setting numpy function
        ## if function should return one constant value
        if self._same_low_high():
            self._samp_n = self._same_low_high()
            self._samp_dist = lambda size: [self._samp_n] * size
            return 0
            
        ## if numpy function
        try:
            setattr(self, '_samp_dist', getattr(np.random, x))            
        except AttributeError:
            raise AttributeError('Distribution "{}" not supported\n'.format(x))

        self._samp_dist = partial(self._samp_dist, **self.samp_dist_params)
        
        try:
            self._samp_dist(size=1)
        except TypeError:
            params = ','.join([':'.join(y) for y in self.samp_dist_params.items()])
            raise TypeError('Params "{}" do not work with distribution "{}"\n'\
                             .format(x, params))        
        return 0
        
        
    
import os, sys
import time
from functools import partial
from collections import Counter
from pprint import pprint
## 3rd party
import numpy as np
import pandas as pd
import dill
from pathos.multiprocessing import ProcessingPool
## amplication
from CommTable import CommTable
from FracTable import FracTable
import BD_Shift
import Utils
from Utils import _table


def binNum2ID(frag_BD_bins, libFracBins):
    """Convert Counter(np.digitize()) for frag_BD  to the fraction BD-min-max.
    Args:
    frag_BD_bins -- dict of counts from np.digitize() binning
    libFracBins -- ordered set of fraction start-ends
    Return:
    dict -- {fractionID : fragment_counts}
    """    
    msg = '{0:.3f}-{1:.3f}'
    n_bins = len(libFracBins)
    out = {}
    for k,v in frag_BD_bins.items():        
        if k < 1:
            binStart = -np.inf
            binEnd = libFracBins[0]
        elif k > len(libFracBins) -1:
            binStart = libFracBins[-1]
            binEnd = np.inf
        else:
            binStart = libFracBins[k-1]
            binEnd = libFracBins[k]
        binRange = msg.format(binStart, binEnd) 
        out[binRange] = str(v)

    return out


def _all_empty_bins(libFracBins):
    """All possible fractions bins are set to zero
    abundance for taxon.
    Args:
    libFracBins -- ordered set of fraction start-ends
    Returns:
    dict -- {bin_ID:'0'}
    """
    n_bins = len(libFracBins)
    return {x:'0' for x in xrange(n_bins+1)}


def _sample_BD_kde(BD_KDE, libID, taxon_name, size):
    """Sampling from buoyant density KDE
    Args:
    BD_KDE -- {library:{taxon:scipy_KDE}}
    libID -- library ID
    taxon_name -- taxon name 
    size -- how many to sample
    Return:
    array -- [fragment BD values]
    """
    try:
        BD_KDE[libID]
    except KeyError:        
        msg = 'Cannot find library->"{}"'
        raise KeyError, msg.format(libID)
    try:
        frag_BD = BD_KDE[libID][taxon_name].resample(size=size)[0]
    except KeyError:
        msg = 'Cannot find taxon->"{}" in lib->"{}"'
        raise KeyError, msg.format(taxon_name, libID)
    except AttributeError:
        frag_BD = None #np.array([np.NaN])
    return frag_BD


def _bin_BD(BD_KDE, libFracBins, taxonAbsAbund, libID, taxon_name,
            maxsize=10000000):
    """Sampling BD KDE and binning values into fractions.
    Parsing binning into chunks to prevent memory errors.
    Args:
    BD_KDE -- {library:{taxon:scipy_KDE}}
    libFracBins -- list of gradient bins
    libID -- library ID
    taxon_name -- taxon name
    maxsize -- max number of BD values to bin at once.
    Return:
    Counter object -- {fraction_ID : count}
    """

    abund_remain = taxonAbsAbund
    frag_BD_bins = Counter()
    while abund_remain > 0:
        if abund_remain < maxsize:
            maxsize = abund_remain
        frag_BD_bins.update(Counter(np.digitize(
            _sample_BD_kde(BD_KDE, libID, taxon_name, maxsize), 
            libFracBins)))        
        abund_remain -= maxsize
        
    return frag_BD_bins


def sim_OTU(x, comm_tbl, libID, BD_KDE, libFracBins, maxsize):
    """Simulate OTU by sampling BD KDE and binning
    values into gradient fractions.
    Args:
    x -- [taxon_idx,taxon_name]
      taxon_idx -- taxon number
      taxon_name -- taxon name 
    comm_tbl -- comm table object
    libID -- library ID
    BD_KDE -- {library:{taxon:scipy_KDE}}
    libFracBins -- list of gradient bins
    maxsize -- max number of BD values to bin at once.    
    Return:
    array -- [fragment BD values]
    """
    taxon_idx,taxon_name = x
    sys.stderr.write('  Processing taxon: "{}"\n'.format(taxon_name))
    

    # taxon abundance based on comm file
    taxonAbsAbund = comm_tbl.get_taxonAbund(libID=libID, 
                                            taxon_name=taxon_name,
                                            abs_abund=True)
    taxonAbsAbund = int(round(taxonAbsAbund[0], 0))
    sys.stderr.write('   taxon abs-abundance:  {}\n'.format(taxonAbsAbund))


    # sampling from BD KDE
    if taxonAbsAbund == 0:
        # zero counts for all bins
       frag_BD_bins =  _all_empty_bins(libFracBins)
    elif taxonAbsAbund < 0:
        raise ValueError, 'Taxon abundance cannot be negative'
    else:
        try:
            frag_BD_bins = _bin_BD(BD_KDE, libFracBins, taxonAbsAbund,
                                   libID, taxon_name, maxsize=maxsize)
        except ValueError:
            # zero counts for all bins
            msg = '   WARNING: No bins for taxon; likely caused by no BD KDE\n'
            sys.stderr.write(msg)
            frag_BD_bins =  _all_empty_bins(libFracBins)

    frag_BD_bins = binNum2ID(frag_BD_bins, libFracBins)

    return [taxon_name, frag_BD_bins] 


def _get_BD_range(x):
    """Getting the BD range from a fraction ID 
    (eg., "1.710-1.716").
    Args:
    x -- fraction ID string
    Returns:
    tuple -- BD start, middle end
    """

    if x.startswith('-'):
        [start,_,end] = x.rpartition('-')
    else:
        [start,_,end] = x.partition('-')            
        
    if start.startswith('-inf'):
        end = round(float(end),3)
        mid = end
    elif end.endswith('inf'):
        start = round(float(start),3)
        mid = start
    else:
        start = round(float(start),3)
        end = round(float(end),3)
        mid = round((end - start)/2 + start,3)
        
    return start, mid, end
    

    
def main(uargs):
    """Main function for making OTU table.
    """
    # args formatting 
    try:
        uargs['--abs'] = int(float(uargs['--abs']))
    except TypeError:
        msg = '"{}" must be float-like'
        raise TypeError(msg.format(uargs['--abs']))

    # logging
    status = Utils.Status(uargs['--quiet'])
    
    # loading files
    sys.stderr.write('Loading files...\n')
    ## BD kde 
    BD_KDE = Utils.load_kde(uargs['<BD_KDE>'])    
    ## community file
    comm_tbl = CommTable.from_csv(uargs['<communities>'], sep='\t')
    comm_tbl.abs_abund = uargs['--abs']
    ## fraction file
    frac_tbl = FracTable.from_csv(uargs['<fractions>'], sep='\t')

    
    # iter by library:
    sys.stderr.write('Simulating OTUs...\n')
    u_taxon_names = comm_tbl.get_unique_taxon_names()
    OTU_counts = []  # list of all library-specific OTU_count dataframes
    for libID in comm_tbl.iter_libraries():
        sys.stderr.write('Processing library: "{}"\n'.format(libID))            

        # fraction bin list for library
        libFracBins = [x for x in frac_tbl.BD_bins(libID)]
        
        # iter of taxa in parallel
        pfunc = partial(sim_OTU, comm_tbl=comm_tbl, libID=libID, 
                        BD_KDE=BD_KDE, libFracBins=libFracBins, 
                        maxsize=int(uargs['--max']))

        pool = ProcessingPool(nodes=int(uargs['--np']))
        if uargs['--debug']:
            ret = map(pfunc, enumerate(u_taxon_names))
        else:
            ret = pool.amap(pfunc, enumerate(u_taxon_names))
            while not ret.ready():
                time.sleep(2)
            ret = ret.get()        

        # converting to a pandas dataframe
        df = pd.DataFrame([x[1] for x in ret]).fillna(0)
        df['taxon'] = [x[0] for x in ret]        
        df = pd.melt(df, id_vars=['taxon'])
        df.columns = ['taxon','fraction', 'count']
        df['library'] = libID
        x = df['fraction'].apply(_get_BD_range).apply(pd.Series)
        x.columns = ['BD_min','BD_mid','BD_max']
        df = pd.concat([df, x], axis=2)
        df = df[['library','taxon','fraction','BD_min','BD_mid','BD_max','count']]
        df.sort(['taxon', 'fraction'], inplace=True)

        # Adding to dataframe of all libraries
        OTU_counts.append(df)

    # combining library-specific dataframes and writing out long form of table
    pd.concat(OTU_counts, ignore_index=False).to_csv(sys.stdout, sep='\t', index=False)



class OTU_table(_table):
    """Subclass of pandas DataFrame
    """

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
        Returns: 
        list -- [min, mean, median, max]
        """
        counts = self.df.groupby(['library','fraction']).sum()['count']
        return [np.min(counts), np.mean(counts), np.median(counts), np.max(counts)]

        
    def subsample(self, no_replace=False):
        """Subsampling from each community.
        Using numpy.random.choice with taxon abundances as weights
        Args:
        no_replace -- subsample without replacement
        Return:
        pandas DataFrame -- subsampled community
        """
        # assertions
        assert hasattr(self, 'samp_dist'), \
            'samp_dist attribute not found'
        assert hasattr(self, 'samp_dist_params'), \
            'samp_dist_params attribute not found'

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
                        sub_comm.loc[:,'fraction'] = fracID
                        BD_min = round(comm['BD_min'].unique()[0], 3)
                        sub_comm.loc[:,'BD_min'] = BD_min
                        BD_mid = round(comm['BD_mid'].unique()[0], 3)
                        sub_comm.loc[:,'BD_mid'] = BD_mid
                        BD_max = round(comm['BD_max'].unique()[0], 3)
                        sub_comm.loc[:,'BD_max'] = BD_max
                    except ValueError:
                        sub_comm = comm.copy()
                        comm.loc[:,'count'] = 0

                df_sub = pd.concat([df_sub, sub_comm])
                        
        df_sub['count'] = df_sub['count'].astype(int)

        cols = ['library','fraction','taxon','BD_min','BD_mid','BD_max','count']
        sort_cols = ['taxon','fraction','library']
        return df_sub.reindex_axis(cols, axis=1).sort(sort_cols)
                
        
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
        """Returns subset of community dataframe.
        Subset selected from libID and fracID args.
        Args:
        libID -- library ID
        fracID -- fraction ID
        Returns:
        pandas.DataFrame
        """
        return self.df.loc[(self.df['library'] == libID) &
                           (self.df['fraction'] == fracID),:]

            
    # iter
    def iter_fractions(self, libID):    
        """iterating through fractions of a certain library.
        Yields a fraction ID.
        Args:
        libID -- library ID
        """
        df_sub = self.df.loc[self.df['library'] == libID]
        for fracID in df_sub['fraction'].unique():
            yield fracID
        
    
    # properties/setters
    @property 
    def is_long(self):
        """Returns true if the table is in 'long' format
        """
        long_cols = ['library','fraction']
        return all([x in self.df.columns.values for x in long_cols])
    @property
    def samp_dist_params(self):
        return self._samp_dist_params
    @samp_dist_params.setter
    def samp_dist_params(self, x):
        self._samp_dist_params = x


    @property
    def columns(self):
        return self.df.columns
    @columns.setter
    def columns(self, x):
        self.df.columns = x

        
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
        
        
    

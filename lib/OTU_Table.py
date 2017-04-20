import os, sys
import time
import math
from functools import partial
from collections import Counter
from pprint import pprint
## 3rd party
import numpy as np
from numpy.random import choice
import pandas as pd
import dill
from pathos.multiprocessing import ProcessingPool
## amplication
from CommTable import CommTable
from FracTable import FracTable
import BD_Shift
import Utils
from Utils import _table
from Utils import random_walk_var_step


def binNum2ID(frag_BD_bins, libFracBins):
    """Convert Counter(np.digitize()) for frag_BD  to the fraction BD-min-max.

    Parameters
    ----------
    frag_BD_bins : dict 
        counts from np.digitize() binning
    libFracBins : ordered set
        fraction start-ends

    Returns
    -------
    dict : {fractionID : fragment_counts}
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

    Parameters
    ----------
    libFracBins : ordered set
        fraction start-ends

    Returns
    -------
    dict : {bin_ID:'0'}
    """
    n_bins = len(libFracBins)
    return {x:'0' for x in xrange(n_bins+1)}


def _sample_BD_kde(BD_KDE, libID, taxon_name, size):
    """Sample from buoyant density KDE.
    
    Parameters
    ----------
    BD_KDE : scipy KDE object
        KDE object for taxon
    libID : str
        library ID
    taxon_name : str
        taxon name 
    size : int
        sample size
    Return:
    numpy.array : fragment BD values
    """
    try:
        frag_BD = BD_KDE.resample(size=size)[0]
    except AttributeError:
        frag_BD = None
    return frag_BD


def _bin_BD(BD_KDE, libFracBins, taxonAbsAbund, libID, taxon_name,
            maxsize=10000000):
    """Sampling BD KDE and binning values into fractions.
    Parsing binning into chunks to prevent memory errors.
    
    Parameters
    ----------
    BD_KDE : scipy KDE object
    libFracBins : list
        gradient bins
    libID : str
        library ID
    taxon_name : str
        taxon name
    maxsize : int
        max number of BD values to bin at once.

    Returns
    -------
    Counter object : {fraction_ID : count}; empty Counter if KDE is None
    """
    abund_remain = taxonAbsAbund
    frag_BD_bins = Counter()
    while abund_remain > 0:
        if abund_remain < maxsize:
            maxsize = abund_remain
            
        vals = _sample_BD_kde(BD_KDE, libID, taxon_name, maxsize)
        
        if vals is not None:
            frag_BD_bins.update(Counter(np.digitize(vals, libFracBins)))        
        else:
            break

        abund_remain -= maxsize
        vals = None
    return frag_BD_bins


def sim_OTU(x, comm_tbl, libID, libFracBins, maxsize):
    """Simulate OTU by sampling BD KDE and binning
    values into gradient fractions.

    Parameters
    ----------
    x : tuple
        (taxon_idx,taxon_name,KDE)
        * taxon_idx -- taxon number
        * taxon_name -- taxon name 
        * KDE -- scipy KDE for taxon
    comm_tbl : comm table object
    libID : str
        library ID
    libFracBins : list 
        gradient bins
    maxsize : int
        max number of BD values to bin at once.    

    Returns
    -------
    numpy.array : [fragment BD values]
    """
    assert len(x) >=3, 'x must be: (taxon_id, taxon_name, kde)'
    taxon_idx,taxon_name,BD_KDE = x

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

    Parameters
    ----------
    x : str
        fraction ID

    Returns
    -------
    tuple -- BD start, middle end
    """
    if x.startswith('-'):
        [start,_,end] = x.rpartition('-')
    else:
        [start,_,end] = x.partition('-')            
        
    if start.startswith('-inf'):
        end = round(float(end),3) - 0.001
        end = round(end,3)
        mid = end
    elif end.endswith('inf'):
        start = round(float(start),3) + 0.001
        start = round(start, 3)
        mid = start
    else:
        start = round(float(start),3)
        end = round(float(end),3)
        mid = round((end - start)/2 + start,3)
        
    return start, mid, end
    

def _get_KDEs_for_libID(KDEs, KDE_type, libID ):
    """Parse out dict of KDE objects for just libID.
    Parsing depends on the KDE type
    """
    err_msg = 'Cannot find KDEs for library->"{}"'
    if KDE_type == 1:
        KDE = {t:k for t,k in KDEs}
    elif KDE_type == 2:
        KDE = KDEs
    elif KDE_type == 3:
        try: 
            KDE = KDEs[libID]
        except KeyError:
            ## kde library not found, duplicating KDE
            raise KeyError, err_msg.format(libID)            
    elif KDE_type == 4:
        try:
            KDE = Utils.load_kde(KDEs[libID])
        except KeyError:
            ## kde library not found, duplicating KDE
            raise KeyError, err_msg.format(libID)            
    else:
        raise ValueError, 'KDE object type not recognized'
    return KDE

    
def main(uargs):
    """Main function for making OTU table.
    
    Parameters
    ----------
    uargs : dict
        See ``OTU_table`` subcommand.
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
    BD_KDE_all = Utils.load_kde(uargs['<BD_KDE>'])    
    BD_KDE_all_type = Utils.KDE_type(BD_KDE_all)
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
        
        # dict of KDEs for library (libID)
        BD_KDE = _get_KDEs_for_libID(BD_KDE_all, BD_KDE_all_type, libID)

        # fraction bin list for library
        frac_bins = frac_tbl.BD_bins(libID)
        assert len(frac_bins) > 0, 'No fractions for library "{}"'.format(libID)
        libFracBins = [x for x in frac_bins]
        
        # iter of taxa in parallel
        pfunc = partial(sim_OTU, 
                        comm_tbl=comm_tbl, 
                        libID=libID, 
                        libFracBins=libFracBins, 
                        maxsize=int(uargs['--max']))

        pool = ProcessingPool(nodes=int(uargs['--np']))
        if uargs['--debug']:
            ret = map(pfunc,[(i,taxon,BD_KDE[taxon]) for
                             i,taxon in enumerate(u_taxon_names)])
        else:
            ret = pool.amap(pfunc,[(i,taxon,BD_KDE[taxon]) for
                                   i,taxon in enumerate(u_taxon_names)])
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
        df = df[['library','taxon','fraction',
                 'BD_min','BD_mid','BD_max','count']]
        df.sort_values(by=['taxon', 'fraction'], inplace=True)

        # Adding to dataframe of all libraries
        OTU_counts.append(df)

    # combining library-specific dataframes
    df_comb = pd.concat(OTU_counts, ignore_index=False)

    # calculating taxon relative abundances 
    df_comb['count'] = df_comb['count'].astype('int')
    f = lambda x: x / x.sum()
    cols = ['library','fraction']
    df_comb['rel_abund'] = df_comb.groupby(cols).transform(f)['count']

    # writing out long form of table
    df_comb.sort_values(by=['library','taxon','BD_mid'], inplace=True)
    df_comb.to_csv(sys.stdout, sep='\t', index=False)



class OTU_table(_table):
    """Subclass of pandas DataFrame.
    """

    def __init__(self, *args, **kwargs):
        _table.__init__(self, *args, **kwargs)
        self._table_columns = self.df.columns.tolist()

    def set_samp_dist(self, samp_dist, samp_dist_params):
        """Setting subsampling size distribution & params.

        Parameters
        ----------
        samp_dist : str 
            a numpy.random distribution attribute
        samp_dist_params : dict
            dict of params for the distribution function
        """
        self.samp_dist_params = samp_dist_params        
        self.samp_dist = samp_dist


    def get_comm_size_stats(self):
        """Getting stats on the size of each community.

        Returns
        -------
        list : [min, mean, median, max]
        """
        counts = self.df.groupby(['library','fraction']).sum()['count']
        return [np.min(counts), np.mean(counts), 
                np.median(counts), np.max(counts)]


        
    def _frac_size_list(self, libID, max_tries=1000):
        """Get a list of fraction community total abundances (sizes).
        
        Parameters
        ----------
        libID : str
            library ID
        max_tries : int
            max tries to select a value btw min_size & max_size

        Returns
        -------
        list : sample sizes (ints)
        """
        samp_sizes = []
        msg = 'Exceeded tries to select a sample size within the range for {}:{}'
        for fracID in self.iter_fractions(libID=libID):
            tries = 0
            while 1:
                samp_size = self.samp_dist(size=1)[0]
                samp_size = int(round(samp_size, 0))
                if samp_size >= self.min_samp_size and \
                   samp_size <= self.max_samp_size:
                    break
                else:
                    tries += 1
                if tries >= max_tries:
                    raise ValueError, msg.format(libID, fracID)
            samp_sizes.append(samp_size)
        return samp_sizes


    def _samp_probs(self, comm, rel_abund=False):
        """Get subsampling weights from taxon counts or relative abundances.
        If counts used, the relative abundances are calculated and used
        as weights.

        Parameters
        ----------
        comm : comm object
            OTU table of one community 
        rel_abund : bool
            use 'rel_abund' column in table; else use 'count'

        Returns
        -------
        list : relative abundances (weights)
        """
        if rel_abund is True:
            rel_abunds = comm['rel_abund']
        else:
            counts = comm['count']
            if self.base is not None:
                counts = [math.log(x+1,self.base) for x in counts]   
            
            rel_abunds = counts / np.sum(counts)                    
        
        # assertion
        total_rel_abund = round(np.sum(rel_abunds),5)
        msg = 'Probabilities = {}, but should = 1'
        assert total_rel_abund == 1.0, msg.format(total_rel_abund)
        
        return rel_abunds
    
        
    def subsample(self, no_replace=False, walk=0, 
                  min_size=0, max_size=None, base=None):
        """Subsample from each community (fraction).
        Using numpy.random.choice with taxon abundances as weights

        Parameters
        ----------
        no_replace : bool
            subsample without replacement
        walk : int
            order values by a random walk with `walk` as max step size
        min_size, max_size : int
            min|max sample size
        base : float
            log base for transforming taxon counts 
            (used for subsampling probabilities)

        Returns
        -------
        pandas DataFrame : subsampled community table
        """
        # attributes
        self.walk = walk
        self.min_samp_size = min_size
        self.max_samp_size = max_size
        self.base = base

        # assertions
        msg = '{} attribute not found'
        assert hasattr(self, 'samp_dist'), msg.format('samp_dist')
        assert hasattr(self, 'samp_dist_params'), msg.format('samp_dist_params')
        msg = 'min_size must be >= 0'
        assert self.min_samp_size >= 0, msg
        msg = 'max_size must be >= min_size'
        assert self.max_samp_size >= self.min_samp_size, msg

        # counting all taxa
        all_taxa = Counter({x:0 for x in self.iter_taxa()})        

        # subsampling
        #samp_cnt = 0
        df_sub = pd.DataFrame(columns=self.df.columns)
        for libID in self.iter_libraries():
            # making list of fraction sizes 
            samp_sizes = self._frac_size_list(libID)

            # applying autocorrleation via random walk (if needed)
            if walk > 0:
                samp_sizes = random_walk_var_step(samp_sizes, walk)
                
            # assertion on number of fractions
            nfracs = len([x for x in self.iter_fractions(libID=libID)])
            msg = 'Number of sample sizes ({}) != number of fractions ({})'
            assert len(samp_sizes) == nfracs, msg.format(len(samp_sizes),
                                                             nfracs)

            # subsampling; weighted by relative abundances in sample
            for samp_cnt, fracID in enumerate(self.iter_fractions(libID=libID)):
                # single community
                comm = self.get_comm(libID, fracID)
                sub_comm = None

                # if no counts for community
                if np.sum(comm['count']) < 1:
                    sub_comm = comm.copy()
                else:
                    # size to subsample
                    samp_size = samp_sizes[samp_cnt]
                    
                    # sampling probabilities (relative abundances)                    
                    rel_abunds = self._samp_probs(comm, rel_abund=True)

                    # sample from taxa list
                    try:
                        sub_comm = Counter(choice(comm['taxon'],
                                                  size=samp_size,
                                                  replace= not no_replace,
                                                  p=rel_abunds))
                                                
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
        df_sub['rel_abund'] = self._norm_counts(df_sub)

        cols = ['library','fraction','taxon',
                'BD_min','BD_mid','BD_max',
                'count', 'rel_abund']
        sort_cols = ['library','taxon','BD_mid']
        return df_sub.reindex_axis(cols, axis=1).sort_values(by=sort_cols)
                
        
    def _same_low_high(self, ret=False):
        """Check for same 'low' and 'high' parameter values.
        If same (and ret=False), returns True, else False.
        If ret=True, the 'low|high' value is returned.
        """
        try:
            same = self.samp_dist_params['high'] == self.samp_dist_params['low']
            if ret == True:
                return self.samp_dist_params['high']
            else: 
                return same
        except KeyError:
            return False
        except AttributeError:
            return False

    def transform_by_group(self, f, sel_index, val_index=['values'], 
                           groups=['library','fraction'], inplace=True):
        """Use the pandas.DataFrame.groupby.transform() function
        to apply a function to the column(s) selected ('sel_index').
        The function must return the same sized object.
        By default, the function will be applied to each community (sample).

        Parameters
        ----------
        f : function
            Applied function
        sel_index : pandas column index
            column name(s) for apply the function to.
        val_index : pandas column index
            column name(s) for assigning resulting values.
        groups : pandas column index (list)
            Which columns to use for grouping. 
            By community(sample) = ['library', 'fraction']
            By taxon (& library) = ['library', 'taxon']
        inplace : bool
            in-place edit of the OTU table.
            The function must not be aggregating (eg., `sum`)
        """
        try:
            self.df.groupby(groups)
        except KeyError:
            msg = "You must select columns in the OTU table:\n  {}"
            cols = '\n  '.join(self.df.columns)
            sys.exit(msg.format(cols))

        x = self.df.groupby(groups)[sel_index].transform(f)        

        if inplace == True:
            self.df[val_index] = x
            return None
        else:
            return x.reset_index()

            
    def apply_by_group(self, f, val_index='values', 
                       groups=['library','fraction'], 
                       inplace=True):
        """Apply a function to each grouping in OTU table.
        By default, the function will be applied to each community (sample).

        Apply will pass each row of the grouped dataframes to the function.
        To select columns, use something like: `lambda x: x['A'] + x['B']`

        Parameters
        ----------
        f : function
            Applied function
        val_index : string
            column name for assigning resulting values.
        groups : pandas column index (list)
            Which columns to use for grouping. 
            By community(sample) = ['library', 'fraction']
            By taxon (& library) = ['library', 'taxon']
        inplace : bool
            in-place edit of the OTU table.
            The function must not be aggregating (eg., `sum`)
        """
        try:
            self.df.groupby(groups)
        except KeyError:
            msg = "You must select columns in the OTU table:\n  {}"
            cols = '\n  '.join(self.df.columns)
            sys.exit(msg.format(cols))
            
        x = self.df.groupby(groups).apply(f)

        if inplace == True:
            self.df[val_index] = x.tolist()
            return None
        else:
            try:
                x.index.names = list(groups) + [val_index]
                x = x.reset_index()
            except ValueError:
                x = x.reset_index()
                ncol = len(x.columns)
                x.columns = x.columns[:ncol-1].tolist() + [val_index]
            return x
      

    def apply_each_taxon(self, f, val_index):
        """Apply a function to each taxon (all entries) in OTU table.
        In-place edit of OTU table.

        Parameters
        ----------
        f : function
            Applied function
        val_index : pandas column index
            column(s) of resulting values        
        """
        self.df[val_index] = self.df[val_index].apply(f)


    def _norm_counts(self, df, sel_index='count'):
        """Normalize OTU counts by total count for community.

        Parameters
        ----------
        df : pandas dataframe
            OTU table
        sel_index :  column with OTU count values

        Returns
        -------
        pandas series of normalized values
        """
        f = lambda x: np.nan_to_num(x / np.sum(x))
        return df.groupby(['library','fraction'])[sel_index].transform(f)


    def add_rel_abund(self, sel_index='count', val_index='rel_abund'):
        """Adding relative abundances (fractions) to OTU table.
        In-place edit of OTU table: new/editted column = 'rel_abund'

        Parameters
        ----------
        sel_index : pandas column index
            column(s) use for for calculations
        val_index : pandas column index
            column(s) of resulting values        
        """        
        self.df[val_index] = self._norm_counts(self.df, sel_index)    


    def adjust_abs_abund(self, rel_index, abs_index):
        """Adjust the absolute values by relative values.
        The abs_index column is edited in place.

        Parameters
        ----------
        sel_index : pandas column index
            column(s) use for for calculations
        val_index : pandas column index
            column(s) of resulting values        
        """
        abs_sum = np.sum(self.df.loc[:,abs_index])
        x = self.df.loc[:,rel_index] * float(abs_sum)
        f = lambda x : round(x, 0)
        self.df[abs_index] = x.apply(f, axis=1).astype('int')

        
    def rm_columns(self, index):
        """Delete/remove columns from the OTU table.
        Columns removed in place.

        Parameters
        ----------
        f : function
            Applied function
        index : pandas column index
            column(s) to drop
        """
        self.df.drop(index, axis=1, inplace=True)

    def rename_columns(self, *args, **kwargs):
        """Renaming columns of OTU table.

        Parameters
        ----------
        See pandas.DataFrame.rename
        """
        self.df.rename(*args, **kwargs)

    def reset_columns(self):
        """Reseting OTU table column order based on initial ordering.
        """
        self.df = self.df[self._table_columns]


    def get_comm(self, libID, fracID):        
        """Returns subset of community dataframe.
        Subset selected from libID and fracID args.

        Parameters
        ----------
        libID : str
            library ID
        fracID : str
            fraction ID
        Returns
        -------
        pandas.DataFrame : subset of community dataframe
        """
        return self.df.loc[(self.df['library'] == libID) &
                           (self.df['fraction'] == fracID),:]


    def select(self, columns):
        """Get certain columns from the OTU table.
        
        Parameters
        ----------
        columns : pandas column index
        """
        return self.df[columns]

    def sort_values(self, **kwargs):
        self.df.sort_values(**kwargs)


    def to_csv(self, *args, **kwargs):
        """Write OTU table as CSV.

        Parameters
        ----------
        See pandas.to_csv
        """
        self.df.to_csv(*args, **kwargs)

            
    # iter
    def iter_fractions(self, libID):    
        """iterating through fractions of a certain library.
        
        Parameters
        ----------
        libID : str
            library ID

        Yields
        ------
        str : fraction ID
        """
        df_sub = self.df.loc[self.df['library'] == libID]
        for fracID in df_sub['fraction'].unique():
            yield fracID


    def merge(self, *args, **kwargs):
        """pandas `merge` for the OTU table dataframe.
        WARNING: in-place edit.
        """
        self.df = self.df.merge(*args, **kwargs)

    def drop(self, *args, **kwargs):
        """Remove rows/colums; WARNING: in-place edit
        """
        self.df = self.df.drop(*args, **kwargs)
    
    # properties/setters
    @property 
    def is_long(self):
        """Return true if the table is in 'long' format"""
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
            self._samp_n = self._same_low_high(ret=True)
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
            params = ','.join([str(x) + ':' + str(y)  for x,y 
                               in self.samp_dist_params.items()])
            raise TypeError('Params "{}" do not work with distribution "{}"\n'\
                             .format(x, params))        
        return 0
                
    @property
    def min_samp_size(self):
        return self._min_samp_size
    @min_samp_size.setter
    def min_samp_size(self, x):
        if x is None:
            x = 0
        self._min_samp_size = float(x)
    
    @property
    def max_samp_size(self):
        return self._max_samp_size
    @max_samp_size.setter
    def max_samp_size(self, x):
        if x is None:
            x = np.inf
        self._max_samp_size = float(x)
    
    @property
    def base(self):
        return self._base
    @base.setter
    def base(self, x):
        if x is None:
            self._base = None
        else:
            self._base = float(x)

    

import os, sys
from functools import partial
from collections import Counter
## 3rd party
import parmap
import numpy as np
import pandas as pd
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
        out[binRange] = v

    return out

def sample_BD_kde(BD_KDE, libID, taxon_name, size):
    try:
        frag_BD = BD_KDE[libID][taxon_name].resample(size=size)[0]
    except KeyError:
        msg = 'Cannot find lib->"{}":taxon->"{}"'
        raise KeyError, msg.format(libID, taxon_name)
    except AttributeError:
        frag_BD = np.array([np.NaN])
    return frag_BD

    
def main(args):
    # args formatting 
    try:
        args['--abs'] = int(float(args['--abs']))
    except TypeError:
        msg = '"{}" must be float-like'
        raise TypeError(msg.format(args['--abs']))

    # logging
    status = Utils.Status(args['--quiet'])
    
    # loading files
    sys.stderr.write('Loading files...\n')
    ## BD kde 
    BD_KDE = Utils.load_kde(args['<BD_KDE>'])    
    ## community file
    comm_tbl = CommTable.from_csv(args['<communities>'], sep='\t')
    comm_tbl.abs_abund = args['--abs']
    ## fraction file
    frac_tbl = FracTable.from_csv(args['<fractions>'], sep='\t')

    
    # iter by library:
    sys.stderr.write('Simulating OTUs...\n')
    u_taxon_names = comm_tbl.get_unique_taxon_names()
    OTU_counts = []  # list of all library-specific OTU_count dataframes

    for libID in comm_tbl.iter_libraries():
        sys.stderr.write('Processing library: "{}"\n'.format(libID))            

        # fraction bin list for library
        libFracBins = [x for x in frac_tbl.BD_bins(libID)]

        
        # creating a dataframe of fraction bins
        func = lambda x: '{0:.3f}-{1:.3f}'.format(libFracBins[x-1],libFracBins[x])
        fracBins = [func(i) for i in xrange(len(libFracBins))][1:]        
        lib_OTU_counts = pd.DataFrame({'fractions':fracBins})
        
        # iter by taxon:
        for (taxon_idx,taxon_name) in enumerate(u_taxon_names):
            sys.stderr.write('  Processing taxon: "{}"\n'.format(taxon_name))

            # taxon abundance
            taxonAbsAbund = comm_tbl.get_taxonAbund(libID=libID, 
                                                    taxon_name=taxon_name,
                                                    abs_abund=True)
            taxonAbsAbund = int(round(taxonAbsAbund[0], 0))
            sys.stderr.write('   N-fragments:   {}\n'.format(taxonAbsAbund))
            
            # sampling from BD KDE
            if taxonAbsAbund == 0:
                frag_BD_bins = {x:0 for x in lib_OTU_counts['fractions']}                
            elif taxonAbsAbund < 0:
                raise ValueError, 'Taxon abundance cannot be negative'
            else:
                try:
                    frag_BD_bins = Counter(np.digitize(
                        sample_BD_kde(BD_KDE, 
                                      libID, 
                                      taxon_name, 
                                      taxonAbsAbund), 
                        libFracBins))
                except ValueError:
                    raise ValueError, 'No values returned from BD KDE'
                frag_BD_bins = binNum2ID(frag_BD_bins, libFracBins)

                       
            # converting to a pandas dataframe
            df = pd.DataFrame(frag_BD_bins.items())
            df.columns = ['fractions',taxon_name]
            
            # adding to dataframe
            df.iloc[:,1] = df.applymap(str).iloc[:,1]   # must convert values to dtype=object
            lib_OTU_counts = pd.merge(lib_OTU_counts, df, how='outer', on='fractions')


        # formatting completed dataframe of OTU counts for the library
        lib_OTU_counts.fillna(0, inplace=True)
        lib_OTU_counts = pd.melt(lib_OTU_counts, id_vars=['fractions'])
        lib_OTU_counts.columns = ['fractions', 'taxon', 'count']
        lib_OTU_counts['library'] = libID
        lib_OTU_counts = lib_OTU_counts[['library','fractions','taxon','count']]
        lib_OTU_counts.sort(['fractions'])


        # making dict of library-specific dataframes
        OTU_counts.append(lib_OTU_counts)

    # combining library-specific dataframes and writing out long form of table
    pd.concat(OTU_counts, ignore_index=False).to_csv(sys.stdout, sep='\t', index=False)




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
        return df_sub.reindex_axis(['library','fractions','taxon','count'], axis=1)\
            .sort(['taxon','fractions','library'])
                
        
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
        
        
    

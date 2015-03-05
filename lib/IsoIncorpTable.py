# import
## batteries
from collections import defaultdict
## 3rd party
#import pymix.mixture as mixture
import mixture
## applicaton
from Utils import _table


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
                # assert start <= end
                if startParam > endParam:
                    params['start'] = endParam
                    params['end'] = startParam
                elif startParam == endParam:  
                    if startParam >= 100:
                        params['start'] -= 1e-5
                    else:
                        params['end'] += 1e-5

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

        return incorpFunc.sampleSet(n_samples)

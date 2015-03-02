# import
## batteries
import sys, os
import numpy as np
import pandas as pd
import scipy.stats as stats
from StringIO import StringIO
from random import randrange, sample, shuffle
import re
from functools import partial
from itertools import chain
from operator import itemgetter
from collections import OrderedDict

## 3rd party
from configobj import ConfigObj, flatten_errors
from validate import Validator


# utility functions
def str2dict(s):
    """Parsing string (format: 'item:value,item:value')
    to create a dict object
    """
    
    if hasattr(s, 'split'):
        l = re.split('[:,]', s)
        try:
            return {k.lower():float(v) for k,v in zip(l[0::2],l[1::2])}
        except TypeError:
            msg = 'distribution parameter values must be ints or floats.'
            raise TypeError(msg)
    else:
        return s
        
        
def random_insert_seq(l, seq):
    """Insert seq items at random locations in list.
    Args:
    l  -- target list 
    seq -- iterable of items to insert
    Return:
    list of with values randomly inserted
    """
    insert_locs = sample(xrange(len(l) + len(seq)), len(seq))
    inserts = dict(zip(insert_locs, seq))
    inputs = iter(l)
    return [inserts[pos] if pos in inserts else next(inputs)
              for pos in xrange(len(l) + len(seq))]


    
class _Comm(object):
    """Parent class for other classes in the module
    """
    
    def __init__(self):
        pass
    
    # property/setter
    @property
    def abund_dist(self):
        return self._abund_dist
    @abund_dist.setter
    def abund_dist(self, x):
        self._abund_dist = str(x)
    @property
    def richness(self):
        return self._richness
    @richness.setter
    def richness(self, x):
        x = float(x)
        if x <= 1:
            x = len(self.taxon_pool) * x
        self._richness = int(round(x,0))


        
class SimComms(_Comm):
    """Simulating communities"""

    def __init__(self, taxon_list, perm_perc, shared_perc,
                 richness, abund_dist, abund_dist_params,
                 n_comm, config=None,
                 *args, **kwargs):
        _Comm.__init__(self, *args, **kwargs)
        
        self._load_taxon_list(taxon_list)
        self.perm_perc = perm_perc
        self.shared_perc = shared_perc
        self.richness = richness
        self.abund_dist = abund_dist
        self.abund_dist_params = str2dict(abund_dist_params)
        self.config = config
        self.n_comm = n_comm
                    
        
        # loading config; setting community parameters
        if config is not None:
            self.comm_params = self._load_config()
        else:
            self.comm_params = dict()
        self._set_comm_params()

        # shared taxa
        self._set_shared_taxa()
                
            
    def _get_configspec(self, strIO=True):
        """Return configspec set for instance.
        Args:
        strIO -- return configspec as a StringIO instance
        """
        configspec = """
        [__many__]
            richness = float(0,inf, default=None)        
            abund_dist = string(default='exponential,1,0.5')
            start = float(default=None)
            end = float(default=None)
            loc = float(default=None)
            scale = float(default=None)
        """
        
        if strIO == True:
            return StringIO(configspec)
        else:
            return configspec

            
    def _load_config(self):
        assert hasattr(self, 'config'), "No config attribute found."

        configspec = self._get_configspec()
        return ConfigObj(self.config, configspec=configspec)        

        
    def _set_comm_params(self):
        """Setting community-specific params including applying global params.
        """
        # adding to comm params if not enough set by config
        n_config_comms = len(self.comm_params.keys())
        n_diff = self.n_comm - n_config_comms
        
        for i in xrange(n_diff):
            self.comm_params[str(n_config_comms + i + 1)] = dict()
        
        for k,v in self.comm_params.items():
            # checking for params
            if ('richness' not in v.keys() or
                v['richness'] is None):
                v['richness'] = self.richness
            if ('abund_dist' not in v.keys() or
                v['abund_dist'] is None):
                v['abund_dist'] = self.abund_dist
            if ('abund_dist_p' not in v.keys() or
                v['abund_dist_p'] is None):
                v['abund_dist_p'] = self.abund_dist_params
            v['abund_dist_p'] = str2dict(v['abund_dist_p'])

                
    def _set_shared_taxa(self):
        """A list of taxa shared among all communities.
        The taxon list (pool) is reduced to just unshared taxa.
        """
        self.shared_taxa = self._drawFromTaxonPool(self.n_shared)

        
    def _load_taxon_list(self, fileName):
        """Loading taxon list file. Taxa order is randomly shuffled.
        Args:
        fileName -- name of taxon file
        """
        self._taxon_pool = []
        if fileName == '-':
            for l in sys.stdin:
                self._taxon_pool.append(l.rstrip())
        else:
            with open(fileName, 'r') as inF:
                for l in inF.readlines():
                    self._taxon_pool.append(l.rstrip())
        shuffle(self._taxon_pool)

        
    def _drawFromTaxonPool(self, n):
        """Drawing from taxon pool, returning n-taxa;
        those taxa are removed from the pool.
        Args:
        n -- number of taxa to draw
        Return:
        list of taxon names
        """
        assert n <= len(self.taxon_pool), \
            'Cannot draw {} taxa from taxon pool'.format(n)
        taxa = self.taxon_pool[:n]
        self.taxon_pool = self.taxon_pool[n:]
        return taxa

                
    def make_comm(self, comm_id):
        """Making a Comm object.
        Args:
        comm_id -- str; community name from comm_params attrib
        """
        # assertions
        comm_id = str(comm_id)
        if not hasattr(self, 'comms'):
            self.comms = OrderedDict() #dict()
            
        try:
            self.comm_params[comm_id]
        except KeyError:
            raise KeyError('Cannot find community ID "{}" in '\
                           'community params\n'.format(comm_id))

        # init comm objects
        self.comms[comm_id] = Comm(comm_id, self)

        
    def write_comm_table(self, Long=True):
        """Joining comm objects into 1 dataframe and printing.
        Args:
        Long -- write table in long format
        """
        df =  pd.concat([x.taxa for x in self.values()],
                        axis=1)
        write_index = True
        df.columns = self.keys()
        if Long == True:
            write_index = False
            # melting
            val_vars = list(df.columns)
            df['taxon'] = df.index
            df = pd.melt(df, id_vars=['taxon'], value_vars=val_vars)
            # ordering columns
            df.columns = ['taxon_name', 'library', 'rel_abund_perc']            
            df = df[['library','taxon_name','rel_abund_perc']]
            # getting rank by community (grouping by community)
            ## TODO
            
        # writing dataframe
        df.to_csv(sys.stdout, sep='\t', na_rep=0,
                  float_format='%.9f', index=write_index)
        
        
    @staticmethod
    def permute(comm, perm_perc):
        """Permute a certain percentage of the taxa abundances.
        Permuting just the indices of the series objects.
        """
        # assertions
        perm_perc = float(perm_perc)
        assert (perm_perc >= 0 and perm_perc <= 100),\
            'perm_perc is not in range [0,100]'
        assert hasattr(comm, 'taxa'), \
            'No "taxa" attribute for comm {}'.format(comm.comm_id)

        # variables
        n_perm = int(round(perm_perc / 100 * comm.n_taxa,0))
        
        # permuting index of comm
        perm_idx = sample(range(comm.n_taxa), n_perm)
        perm_ig = itemgetter(perm_idx)
        n_perm_idx = set(range(comm.n_taxa)) - set(perm_idx)
        if len(n_perm_idx) > 0:
            n_perm_ig = itemgetter(*n_perm_idx)        
            # altering pandas series of taxa & abundances
            comm.taxa.index = random_insert_seq(n_perm_ig(comm.taxa.index),
                                                perm_ig(comm.taxa.index))
        else:
            # altering pandas series of taxa & abundances
            comm.taxa.index = random_insert_seq([],
                                                perm_ig(comm.taxa.index))
            
        
        
    def items(self):
        return self.comms.items()
    def keys(self):
        try:
            return self.comms.keys()
        except AttributeError:
            return np.sort(self.comm_params.keys())
    def values(self):
        return self.comms.values()
        
    # property/setter
    @property
    def n_comm(self):        
        return self._n_comm
    @n_comm.setter
    def n_comm(self, x):
        try:
            self._n_comm = int(x)
        except ValueError:
            raise ValueError('n_comm must be an integer')
    @property
    def perm_perc(self):
        return self._perm_perc
    @perm_perc.setter
    def perm_perc(self, x):
        x = float(x)
        assert (x >= 0 and x <= 100), 'shared_perc must be in range 0-100'
        self._perm_perc = x        
    @property
    def shared_perc(self):
        return self._shared_perc
    @shared_perc.setter
    def shared_perc(self, x):
        x = float(x)
        assert (x >= 0 and x <= 100), 'shared_perc must be in range 0-100'
        self._shared_perc = x
    @property
    def taxon_pool(self):
        return self._taxon_pool
    @taxon_pool.setter    
    def taxon_pool(self, x):
        self._taxon_pool = x

    @property
    def min_richness(self):
        """The minimum richness of any community as defined by comm_params.
        """
        if not hasattr(self, '_min_richness'):
            setattr(self, '_min_richness', None)
        if self._min_richness is None:
            richness_vals = []
            for k,v in self.comm_params.items():
                try:
                    richness_vals.append(int(v['richness']))
                except KeyError:
                    raise KeyError('Cannot find richness attribute for '\
                                   'comm_id "{}"'.format(k))
            self._min_richness = min(richness_vals)
        return self._min_richness
        
                    
    @property
    def n_shared(self):
        """The number of taxa that should be shared;
        defined by shared_perc * min richness of any community.
        """
        if not hasattr(self, '_n_shared'):
            setattr(self, '_n_shared', None)
        if self._n_shared is None:
            self._n_shared = self.min_richness * (self.shared_perc / 100.0)
            self._n_shared = int(round(self._n_shared,0))
        return self._n_shared

    @property
    def n_taxa_remaining(self):
        """The number of taxa that remain in taxa pool.
        """
        if not hasattr(self, '_n_taxa_remaining'):
            setattr(self, '_n_taxa_remaining', None)
        return len(self.taxon_pool)
            
        


class Comm(_Comm):
    """Community class"""

    def __init__(self, comm_id, GradientComms, *args, **kwargs): 
        """
        Args:
        taxon_list -- list of taxa
        richness -- number of taxa in community
        abund_dist -- taxon abundance distribution
        """
        _Comm.__init__(self, *args, **kwargs)
        self.comm_id = comm_id                                                     
        self.params = GradientComms.comm_params[comm_id]
        self.n_shared = GradientComms.n_shared
        self.richness = self.params['richness']
                                                       
        # assertions
        if self.richness > self.n_shared + GradientComms.n_taxa_remaining:
            sys.exit('ERROR: Comm_ID {}\n'\
            '  Community richness is set too high! It is > taxon pool.\n'\
            '  There is not enough taxa to make the desired communities.\n' \
            '  You must reduce richness or increase perc_shared.\n'\
            '  NOTE: shared_perc is based on the community with the min. richness.\n'\
            .format(comm_id))
        
        # selecting additional taxa beyond those shared by all comms
        ## unique taxa inserted randomly in list while perserving shared taxa rank-abund
        n_unique = self.richness - GradientComms.n_shared
        assert n_unique >= 0, 'ERROR: Comm_ID {}: the number ' \
            'of unique taxa is < 0'.format(comm_id)
        self.taxa = random_insert_seq(GradientComms.shared_taxa,
                                      GradientComms._drawFromTaxonPool(n_unique))

        
        # drawing relative abundances from the user-defined distribution
        abund_dist = self._get_abund_dist(self.params['abund_dist'],
                                          self.params['abund_dist_p'])        
        rel_abunds = abund_dist(size=self.n_taxa)
        rel_abunds = np.sort(rel_abunds / sum(rel_abunds) * 100)[::-1]
    
        # making a series for the taxa
        self.taxa = pd.Series(rel_abunds, index=self.taxa)        


    def __repr__(self):
        return self.taxa.__repr__()

        
    def _get_abund_dist(self, dist, params):
        try:
            distFunc = getattr(np.random, dist)
        except AttributeError:
            msg = 'Distribution "{}" is not supported'.format(dist)
            
        try:
            return partial(distFunc, **params)
        except TypeError:
            param_str = [':'.join([str(k),str(v)]) for k,v in params.items()]
            param_str = ','.join(param_str)
            msg = 'Params "{}" do not work with function "{}"'\
                  .format(param_str, dist)


    @property
    def n_taxa(self):
        return len(self.taxa)

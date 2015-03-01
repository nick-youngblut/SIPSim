# import
## batteries
import sys, os
import numpy as np
import pandas as pd
import scipy.stats as stats
from StringIO import StringIO
import random
import re
from functools import partial
from itertools import chain

## 3rd party
from configobj import ConfigObj, flatten_errors
from validate import Validator


#-- workflow --#

## init:
### load args
### load list of taxa
### config file:
#### global: 
##### permute percent -- % of rank-abundances permuted in each library
##### shared percent -- % of taxa shared between libraries
##### richness -- number of taxa in each library (0-1, min to max possible)
##### abund distribution -- 
#### per-library:
##### richness -- number of taxa in the library
##### abund distribution 



#-- community diversity --#
# shared perc:
## based on richness (% of community with smallest richness)
## must have enough taxa in total pool for all unshared taxa
### total unshared needed = sum(richness_comX - shared) 
## randomly select shared taxa from taxa pool
## taxa will have same relative rank-abundance ordering

# community total richness:
## for each Comm: randomly select (without replacemnt) other, non-shared taxa from pool
## taxa randomly inserted into current rank-abundance ordering 

# permuted percent
## for each Comm
## permute rank-abundances of X% of taxa

# abundance distribution (numpy / scipy)
## for defined richness, draw from distribution to get abundance for each rank



#-- OO design --#
# class GradientComms
## global params
## list of taxa
## (opt) list of taxon abundances
## creates Comm class objects based on config file
## for each Comm, either applies global param or sets comm-specific param
## functions for permuting communties 

# class Comm
## community of taxa and relative abundances

def str2dict(s):
    """Parsing string (format: 'item:value,item:value')
    to create a dict object
    """
    if s is None or type(s) is dict:
        return s
    elif type(s) is str:
        l = re.split('[:,]', s)
        try:
            return {k.lower():float(v) for k,v in zip(l[0::2],l[1::2])}
        except TypeError:
            msg = 'distribution parameter values must be ints or floats.'
            raise TypeError(msg)
    else:
        raise TypeError


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
    # @property
    # def abund_dist_params(self):
    #     return self._abund_dist_params
    # @abund_dist_params.setter
    # def abund_dist_params(self, x):
    #     try:
    #         x.items()
    #     except AttributeError:
    #         msg = 'params must be a dict-like object'
    #         raise AttributeError(msg)            
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
        for i in xrange(1, self.n_comm + 1):
            i = str(i)
            # setting community
            try:
                self.comm_params[i]
            except:
                self.comm_params[i] = dict()
            # checking for params
            if ('richness' not in self.comm_params[i].keys() or
                self.comm_params[i]['richness'] is None):
                self.comm_params[i]['richness'] = self.richness
            if ('abund_dist' not in self.comm_params[i].keys() or
                self.comm_params[i]['abund_dist'] is None):
                self.comm_params[i]['abund_dist'] = self.abund_dist
            if ('abund_dist_p' not in self.comm_params[i].keys() or
                self.comm_params[i]['abund_dist_p'] is None):
                self.comm_params[i]['abund_dist_p'] = self.abund_dist
            self.comm_params[i]['abund_dist_p'] = \
                        str2dict(self.comm_params[i]['abund_dist_p'])
            
                
                
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
        random.shuffle(self._taxon_pool)

        
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
            self.comms = dict()
            
        try:
            self.comm_params[comm_id]
        except KeyError:
            raise KeyError('Cannot find community ID "{}" in '\
                           'community params\n'.format(comm_id))

        # init comm objects
        self.comms[comm_id] = Comm(comm_id, self)

        
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
        self.params = GradientComms.comm_params[comm_id]
        self.n_shared = GradientComms.n_shared
        self.richness = self.params['richness']
                                                       
        # assertions
        assert self.richness <= self.n_shared + GradientComms.n_taxa_remaining,\
            'Comm_ID {}: richness is > taxon pool.\n'\
            '  There is not enough taxa to make the desired communities.\n' \
            '  You must reduce community richness or increase perc_shared.'
        
        # selecting additional taxa beyond those shared by all comms
        ## taxa inserted randomly in list
        ## TODO : finish
        n_unique = self.richness - GradientComms.n_shared
        assert n_unique >= 0, 'ERROR: Comm_ID {}: the number ' \
            'of unique taxa is < 0'.format(comm_id)
        #self.taxa = GradientComms.shared_taxa + \
        #            GradientComms._drawFromTaxonPool(n_unique)
        #self.taxa = [0] * richness
        

        # drawing relative abundances from the user-defined distribution
        abund_dist = self._get_abund_dist(self.params['abund_dist'],
                                          self.params['abund_dist_p'])        
        rel_abunds = abund_dist(size=self.n_taxa)
        rel_abunds = rel_abunds / sum(rel_abunds)

        
        # making a series for the taxa
        self.taxa = pd.Series(rel_abunds, index=self.taxa)

        # permuting the rank-abundance of taxa
        ## TODO: finish
        

        
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

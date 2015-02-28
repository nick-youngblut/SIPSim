# import
## batteries
import sys, os
import numpy as np
import pandas as pd
import scipy.stats as stats
from StringIO import StringIO

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



class SimComms(object):
    """Simulating communities"""

    def __init__(self, taxon_list, perm_perc, shared_perc,
                 richness, abund_dist, n_comm, config=None):
        self._perm_perc = perm_perc
        self._shared_perc = shared_perc
        self._richness = richness
        self._abund_dist = abund_dist
        self._config = config
        try:
            self._n_comm = int(n_comm)
        except ValueError:
            raise ValueError('n_comm must be an integer')

        # loading taxon list
        self._load_taxon_list(taxon_list)
        
        # loading config; setting community parameters
        if config is not None:
            self.comm_params = self._load_config()
        else:
            self.comm_params = dict()
        self._set_comm_params()
            
        
        # creating Comm objects
        for i in xrange(self.n_comm):
            self.make_comm(i)

            
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
                

    def _load_taxon_list(self, fileName):
        """Loading taxon list file.
        Args:
        fileName -- name of taxon file
        """
        self._taxon_list = []
        if fileName == '-':
            for l in sys.stdin:
                self._taxon_list.append(l.rstrip())
        else:
            with open(fileName, 'r') as inF:
                for l in inF.readlines():
                    self._taxon_list.append(l.rstrip())
                
        
                
    def make_comm(i):
        """Making a Comm object.
        Args:
        i -- str; community name from comm_params attrib
        """
        # assertions
        i = str(i)
        if not hasattr(self, 'comms'):
            self.comms = []
            
        try:
            self.comm_params[i]
        except KeyError:
            raise KeyError('Cannot find "{}" in community params\n'.format(i))

        # making comm object
        comm = Comm(taxon_list = self.taxon_list,
                    richness = self.comm_params[i]['richness'],
                    abund_dist = self.comm_params[i]['abund_dist'])                    
                    
                

    # properties/setters
    @property
    def n_comm(self):
        return self._n_comm
    @property
    def config(self):
        return self._config
    @property
    def perm_perc(self):
        return self._perm_perc
    @property
    def richness(self):
        return self._richness
    @property
    def abund_dist(self):
        return self._abund_dist
    @property
    def taxon_list(self):
        return self._taxon_list

        
            
class Comm(object):
    """Community class"""

    def __init__(self, richness, abund_dist, taxon_list):
        self._richness = richness
        self._abund_dist = abund_dist

        
        
"""Error distribution functions"""

# import
import sys
import random
## batteries
from functools import partial
## 3rd party
import numpy as np
import pandas as pd


class error_dist(object):
    
    def __init__(self, dist, params):
        """
        Parameters
        ----------
        dist : str
            Name of distribution
        params : dict
            Parameter values {param : value}
        """
        self.params = params
        self.dist = dist


    def sample(self, size, *params, **kparams):        
        """Sample from the error distribution.

        Parameters
        ----------
        size : int
            number of values to sample
        *params, **kparams : args, kargs
            passed to distribution function

        Returns
        -------
        dist_vals : list
            Values drawn from the distribution function.
        """
        try:
            x = self._dist(size=size, *params, **kparams)
            #print [int(y).__class__ for y in x]; sys.exit()
            return [int(y) for y in x]
        except TypeError:
            params = self.params.update(**kparams)
            params = ','.join([str(x) + ':' + str(y)  for x,y 
                               in params.items()])
            msg = 'Params "{}" do not work with distribution "{}"\n'
            raise TypeError(msg.format(params, self._dist_name))        
        

    def _same_low_high(self, ret=False):
        """Check for same 'low' and 'high' parameter values.
        If same (and ret=False), returns True, else False.
        If ret=True, the 'low|high' value is returned.
        """
        try:
            same = self.params['high'] == self.params['low']
            if ret == True:
                return self.params['high']
            else: 
                return same
        except KeyError:
            return False
        except AttributeError:
            return False

    @property
    def dist(self):
        return self._dist
    @dist.setter
    def dist(self, x):
        self._dist_name = x
        # setting numpy function
        ## if function should return one constant value
        if self._same_low_high():
            self._samp_n = self._same_low_high(ret=True)
            self._dist = lambda size: [self._samp_n] * size
        else:
            ## set function
            try:
                setattr(self, '_dist', getattr(np.random, x))            
            except AttributeError:
                msg = 'Distribution "{}" not supported\n'
                raise AttributeError(msg.format(x))
            self._dist = partial(self._dist, **self.params)

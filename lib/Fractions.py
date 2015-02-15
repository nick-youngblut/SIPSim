#import
## batteries
import os,sys
## 3rd party
#import pymix.mixture as mixture
import mixture
import numpy as np


class Fractions(object):
    """Simulated gradient fractions based on theoretical min-max of BD incorporation.
    Fraction size distribution is normal with user-defined params.
    """

    def __init__(self, distribution, params, BD_min, BD_max, frac_min):
        """
        Args:
        distribution -- str; name of distribution for selecting fraction sizes
        params -- dict-like; distribution params passed to distribution function
                  (pymix.mixture distributions)
        BD_min -- Min BD for gradient fractions.
        BD_max -- Max BD for gradient fractions.
        frac_min -- Min size of any fraction.
        """
        BD_min, BD_max, frac_min = map(float, (BD_min, BD_max, frac_min))
        assert BD_min > 0 and BD_max > 0, "BD_min/max must be positive numbers"
        assert BD_max > BD_min, "BD_max must be > BD_min"        
        assert frac_min > 0, "frac_min must be > 0"
        
        self.distribution = distribution
        self.params = params
        self.BD_min = BD_min
        self.BD_max = BD_max
        self.BD_range = BD_max - BD_min
        self.frac_min = frac_min
        
        self._set_distributionFunction()
        
        
    def _set_distributionFunction(self):
        """Setting user-defined distribution (pymix distributions)"""
        psblFuncs = {'uniform' : mixture.UniformDistribution,
                     'normal' : mixture.NormalDistribution}

        try:
            func = psblFuncs[self.distribution.lower()]
        except KeyError:
            msg = 'Distribution "{}" not supported'
            raise KeyError(msg.format(self.distribution))

        try:
            self.distFunc = func(**self.params)
        except TypeError:
            msg = 'Distribution "{}" does not accept all provided params: {}'
            raise TypeError(msg.format(self.distribution, ','.join(self.params.keys())))


    def simFractions(self, libID, max_tries=1000):
        """Simulating the gradient fractions for a library
        Args:
        libID -- library ID (a.k.a. isopycnic gradient ID)
        max_tries -- max number of tries to get BD values in BD_range
        """
        
        # sample & sum to get up to just past (BD_range - BD_max)
        # use truncated last value if value > frac_min
        BD_sums = 0
        BD_vals = []
        tries = 0
        while True:
            fracSize = self.sample_frac_size(max_tries=max_tries)
                    
            BD_sums += fracSize
            BD_vals.append(fracSize)

            BD_prog = self.BD_range - BD_sums
            
            if tries > max_tries:
                msg = 'Exceeded {} tries to make BD fractions. Giving up!\n'
                sys.stderr.write( msg.format(max_tries))
                sys.exit(1)
            elif BD_sums > self.BD_range:
                remains = BD_sums - self.BD_range
                psblLastFrac = round(BD_vals[-1:][0] - remains, 3)

                if psblLastFrac >= self.frac_min:
                    # last fraction set (not exceeding BD-range)
                    BD_vals = BD_vals[:len(BD_vals)-1] + [psblLastFrac]
                    break
                else:
                    # start over
                    BD_sums = 0
                    BD_vals = []
                    tries += 1
                    continue
            else:
                continue

        # making BD-min-max for each fraction
        ## cumulative addition
        BD_vals_cumsum = np.array(BD_vals).cumsum()
        BD_max_vals = BD_vals_cumsum + self.BD_min
        BD_min_vals = np.concatenate( (np.array([self.BD_min]), BD_max_vals[:len(BD_max_vals)-1]) )

        return zip(BD_min_vals, BD_max_vals, BD_vals)


    def sample_frac_size(self, max_tries=1000):
        """Sampling from user-defined distribution. The value must be >= frac_min
        """
        tries = 0
        while True:
            val =  round(self.distFunc.sample(), 3)
            if val >= self.frac_min:
                return val
            else:
                tries += 1

            if tries > max_tries:
                msg = 'Exceeded {} tries to make BD fraction from user-defined distribution. Giving up!'
                print msg.format(max_tries)
                sys.exit(1)

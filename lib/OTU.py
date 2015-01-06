# import
## batteries
import functools
import itertools
## 3rd party
import numpy as np
import pandas as pd
import scipy.stats as stats
## application
from Utils import _table


def add_diffusion(frag_gc_len, loc=0,  max_tries=100):
    """Sampling the diffusion distribution function and adding value
    to the provided GC value.
    Asserts that G+C must be between 0 & 100.
    
    Diffusion equation from: Clay 0, Douady CJ, Carels N, Hughes S, Bucciarelli G, Bernardi G (2003.
    Using analytical ultracentrifugation to study compositional variation in vertebrate genomes.
    Eur Biophys J 32:418-426
    
    Args:
    frag_gc_len -- array: [fragment_gc, fragment_len]
    #gc_val -- G+C value to add noise to
    #frag_len -- fragment length (fragment associated with provided G+C value)
    loc -- loc param of the scipy.stats distribution function
    max_tries -- max number of tries to get value in [gn_range_min,gn_range_max]
    Return:
    generator -- float; new GC value
    """
    loc=0
    max_tries=100
    tries = 0
    while True:
        new_gc = frag_gc_len[0] + np.random.normal(size=1, loc=loc, scale=44500/frag_gc_len[1])
        if new_gc >= 0 and new_gc <= 100:
            return new_gc
            
        tries += 1
        if tries >= max_tries:
            print "ERROR: exceeded max tries ({}) to get"
            "G+C value (with diffusion) between 0 & 100".foramt(max_tries)
            sys.exit(1)
        


class OTU_table(object):
    """Class for the simulated OTU table"""
    
    def __init__(self, frac, g_noise, gn_scale, abund_weight, isotope):
        """
        Args:
        frac -- Fractions class instance
        g_noise -- str; gradient noise distribution name (from scipy.stats)
        gn_scale -- float; scale parameter for gradient noise distribution
        abund_weight -- float; isotope incorp abundance-weighting factor
        isotope -- str; name of isotope
        """        
        self.g_noise = g_noise
        self.gn_scale = float(gn_scale)
        self.abund_weight = abund_weight
        self.isotope = isotope.upper()        
        
        
        # setting gradient 'noise' function
        self.g_noise_func = self._set_g_noise_func()
        
        # set isotope theoretical max BD
        self.isotopeMaxBD = self._set_isotopeMaxBD(self.isotope)
        
    
    def _set_g_noise_func(self):
        """Setting the gradient noise function as scipy distribution function.
        Args:
        #dist_name -- name of scipy distribution function
        #scale -- scale parameter of distribution
        """        
        psblFuncs = {'cauchy' : stats.cauchy,
                     'normal' : stats.norm,
                     'uniform' : stats.uniform}

        try:
            func = psblFuncs[self.g_noise]
        except KeyError:
            raise KeyError('Distribution "{}" not supported.'.format(self.g_noise))
        
        return functools.partial(func, scale=self.gn_scale)            
        
        
    def _set_isotopeMaxBD(self, isotope):
        """Setting the theoretical max BD shift of an isotope (if 100% incorporation).
        Args:
        isotope -- str; name of isotope
        """
        psblIsotopes = {'13C' : 0.036,
                        '15N' : 0.016}
        try:
            return psblIsotopes[isotope]
        except KeyError:
            raise KeyError('Isotope "{}" not supported.'.format(isotope))

            

    def sample_g_noise_func(self, gc_val, loc=0, n_samples=1, max_tries=100):
        """Sampling the gradient noise distribution function and adding
        value to the provided GC value
        Asserts that G+C must be between 0 & 100.
        Args:
        gc_val -- GC value to add noise to
        loc -- loc param of the scipy.stats distribution function
        n_samples -- number of samples to draw from distrubtion
        max_tries -- max number of tries to get value in [gn_range_min,gn_range_max]
        Return:
        float -- new GC value
        """
        tries = 0
        while True:
            new_gc = gc_val + self.g_noise_func(loc=loc).rvs(n_samples)[0]
            if new_gc >= 0 and new_gc <= 100:
                return new_gc
                
            tries += 1
            if tries >= max_tries:
                print "ERROR: exceeded max tries ({}) to get"
                "G+C value (with noise) between 0 & 100".foramt(max_tries)
                sys.exit(1)

            
    def checkLibOverlap(self, libList):
        """Checking that libraries in dataframes fully overlap
        Args:
        libList -- list of list of unique libs (each lib list is compared to the others)
        """        
        return len(set(itertools.chain(*libList))) == len(set(libList[0])) 
        
        
    def make_emptyCountTable_OLD(self, taxon_names, fractionIDs):
        """Making a pandas dataframe (taxon x fraction) of zeros.
        Args:
        taxon_names -- iterable; all taxon names
        fractionIDs -- iterable; all fraction IDs
        """
        shape = (len(taxon_names),len(fractionIDs))
        return pd.DataFrame(np.zeros(shape), columns=fractionIDs, index=taxon_names)
        
        
    def get_isotopeMaxBD(self):
        return self.isotopeMaxBD

    @property
    def gn_scale_nonzero(self):
        return self.gn_scale > 0
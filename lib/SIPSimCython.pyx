from __future__ import division
import sys
import random
import numpy as np
cimport numpy as np
from cython.parallel import prange

import SIPSimCpp

DTYPE = np.float
ctypedef np.float_t DTYPE_t

    
def add_diffusion(np.ndarray[DTYPE_t, ndim=2] arr, 
                  float T=298, float B=1.195e9, 
                  float G=7.87e-10, int M=882):
    """Adding diffusion to fragment buoyant density values. 
    See Clay et al., 2003 for more details.
    Args:
    arr -- 2d-array: [[frag_buoyant_density,],[frag_length,]]
    T -- gradient temperature in Kelvin
    B -- beta coefficient
    G -- G coefficient
    M -- molecular weight per pair base pair of dry cesium DNA
    Returns:
    numpy.array -- [BD values that include diffusion error]
    """
    cdef int n = len(arr[0])
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    cdef double diff_error = 0.0
    for i in xrange(n):
        # fragment length must be > 0
        ## if False: selecting another value from array
        while 1:
            if arr[1,i] <= 0:
                arr[1,i] = arr[1,random.randint(0,n-1)]
            else:
                break

        # error from true BD due to diffusion
        diff_error = SIPSimCpp.calc_diffusion_BD(arr[0,i], 
                                                 arr[1,i], 
                                                 T, B, G, M)            

        # true_BD + diffusion_error_BD
        out[i] = arr[0,i] + diff_error
        
    return out


def GC2BD(np.ndarray[DTYPE_t, ndim=1] arr):
    """Converting G+C (% from 0-100) to buoyant density (BD)
    Equation: (GC / 100) * 0.098 + 1.66
    Args:
    arr -- numpy array of G+C values
    Return:
    numpy.array -- [BD values]
    """
    cdef int n = len(arr)
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    for i in xrange(n):
        out[i] = arr[i] / 100 * 0.098 + 1.66
    return out


def add_incorp(np.ndarray[DTYPE_t, ndim=1] frag_BD, 
               incorp_func, double isotopeMaxBD):
    """Adding isotope incorporation BD-shift values to BD-0%-incorp
    values
    Args:
    frag_BD -- 1d numpy array of fragment BD values
    incorp_func -- function that returns a 1d-list with a float
    isotopeMaxBD -- the max BD possible with the selected isotope 
    Returns:
    float -- buoyant density value
    """    
    cdef int n = len(frag_BD)
    
    cdef int i = 0
    cdef double y = 100.0
    cdef double z = 0
    
    for i in xrange(n):
        z = incorp_func.sample()[0]
        frag_BD[i] += z / y * isotopeMaxBD

    return frag_BD




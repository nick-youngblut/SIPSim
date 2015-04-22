from __future__ import division
import numpy as np
cimport numpy as np

import SIPSimCpp

DTYPE = np.float
ctypedef np.float_t DTYPE_t


    
def add_diffusion_wrapper(np.ndarray[DTYPE_t, ndim=2] arr):
    """Adding diffusion to GC values and calculating BD
    Args:
    arr -- 2d-array: [[GC_value,],[frag_length,]]
    """
    cdef int n = len(arr[0])
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i

    cdef int c1 = 100
    cdef double c2 = 0.098
    cdef double c3 = 1.66
    for i in xrange(n):
        out[i] = SIPSimCpp.add_diffusion(arr[0,i], arr[1,i]) / c1 * c2 + c3
    return out



def add_incorp(frag_BD, incorp_func, double isotopeMaxBD):
    """Adding isotope incorporation BD-shift values to BD-0%-incorp
    values
    Args:
    incorp_func -- function that returns a 1d-list with a float
    isotopeMaxBD -- the max BD possible with the selected isotope 
    """    
    cdef int n = len(frag_BD)
    
    cdef int i = 0
    cdef double y = 100.0
    cdef double z = 0
    
    for i in xrange(n):
        z = incorp_func.sample()[0]
        frag_BD[i] += z / y * isotopeMaxBD




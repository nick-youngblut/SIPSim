from __future__ import division
import numpy as np
cimport numpy as np
from cython.parallel import prange

import SIPSimCpp


DTYPE = np.float
ctypedef np.float_t DTYPE_t


    
def add_diffusion_wrapper(np.ndarray[DTYPE_t, ndim=2] arr):    
    cdef int n = len(arr[0])
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    
    for i in xrange(n):
        out[i] = SIPSimCpp.add_diffusion(arr[0,i], arr[1,i]) / 100.0 * 0.098 + 1.66
    return out


def add_incorp(frag_BD,
               np.ndarray[DTYPE_t, ndim=2] incorp_arr,
               double isotopeMaxBD):
    cdef int n = len(incorp_arr)
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    for i in xrange(n):
        frag_BD[i] += incorp_arr[i,1] / 100.0 * isotopeMaxBD
    return frag_BD
    
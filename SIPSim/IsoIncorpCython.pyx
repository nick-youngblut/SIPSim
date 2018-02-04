from __future__ import division
import sys
import random
import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

    

def GC2BD(np.ndarray[DTYPE_t, ndim=1] arr):
    """Convert G+C (% from 0-100) to buoyant density (BD)
    Equation: (GC / 100) * 0.098 + 1.66

    Parameters
    ----------
    arr : numpy.array 
        G+C values

    Returns
    -------
    numpy.array : BD values
    """
    cdef int n = len(arr)
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    for i in xrange(n):
        out[i] = arr[i] / 100 * 0.098 + 1.66
    return out


def add_incorp(np.ndarray[DTYPE_t, ndim=1] frag_BD, 
               incorp_func, double isotopeMaxBD):
    """Adding isotope incorporation BD-shift values to BD-0%-incorp values.

    Parameters
    ----------
    frag_BD : numpy.array
        1d array of fragment BD values
    incorp_func : function
        A function that returns a 1d-list with a float
    isotopeMaxBD : float
        The max BD possible with the selected isotope 

    Returns
    -------
    float : buoyant density value
    """    
    cdef int n = len(frag_BD)
    
    cdef int i = 0
    cdef double y = 100.0
    cdef double z = 0
    
    for i in xrange(n):
        z = incorp_func.sample()[0]
        frag_BD[i] += z / y * isotopeMaxBD

    return frag_BD





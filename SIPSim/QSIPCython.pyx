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


def calc_density_shift(np.ndarray[DTYPE_t, ndim=1] col_treat, 
                       np.ndarray[DTYPE_t, ndim=1] col_control):
    """Calculating BD shift (treatment - control)
    """
    
    assert len(col_treat) == len(col_control)
    cdef int n = len(col_treat)
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    for i in xrange(n):
        out[i] = col_treat[i] - col_control[i]
    return out

    
def BD2GC(np.ndarray[DTYPE_t, ndim=1] col_BD):
    """buoyant density to G+C (fraction). 
    Ref: Birnie and Rickwood, 1978
    """
    cdef double x = 1.66
    cdef double y = 0.098

    cdef int n = len(col_BD)
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    for i in xrange(n):
        out[i] = (col_BD[i] - x) / y
    return out


def BD2GC(np.ndarray[DTYPE_t, ndim=1] col_BD):
    """buoyant density to G+C (fraction). 
    Ref: Hungate et al., 2015. AEM
    """    
    cdef int x = 1
    cdef double y = 0.083506
    cdef double z = 1.646057

    cdef int n = len(col_BD)
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    for i in xrange(n):
        out[i] = (x / y) * (col_BD[i] - z)
    return out


def GC2M_light(np.ndarray[DTYPE_t, ndim=1] col_GC):
    """G+C to M_light (molecular weight of 'light' DNA).
    Ref: Hungate et al., 2015. AEM
    """
    cdef double x = 0.496
    cdef double y = 307.691

    cdef int n = len(col_GC)
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    for i in xrange(n):
        out[i] = x * col_GC[i] + y
    return out


def M_light2GC(double M_light):
    """The inverse of GC2M_light.
    """
    cdef double x = 307.691
    cdef double y = 0.496
    return (M_light - x) / y

    
def M_light2heavyMax(np.ndarray[DTYPE_t, ndim=1] col_M_light, isotope='13C'):
    """G+C to theoretical molecular weight of fully-labeled DNA (M_HEAVYMAXi).
    Ref: Hungate et al., 2015. AEM.
    
    Parameters
    ----------
    col_M_light : np.array (float)
        Molecular weight of DNA in control gradients (light DNA)
    isotope : str
        13C or 18O isotope?
    """
    cdef double x = -0.4987282
    cdef double y = 9.974564
    cdef double z = 12.07747

    cdef int n = len(col_M_light)
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    for i in xrange(n):
        if isotope == '13C':        
            G = M_light2GC(col_M_light[i])
            out[i] = x * G + y + col_M_light[i]
        elif isotope == '18O':
            out[i] = z + col_M_light[i]
        else:
            msg = '"{}" Isotope not supported; only: "13C" or "18O"'
            raise TypeError, msg.format(isotope)
    return out

def calc_M_lab(np.ndarray[DTYPE_t, ndim=1] col_Z, 
               np.ndarray[DTYPE_t, ndim=1] col_W_light, 
               np.ndarray[DTYPE_t, ndim=1] col_M_light):        
    """Calculate the molecular weight of DNA in labeled treatments.
    Ref: Hungate et al., 2015. AEM.
    
    Parameters
    ----------
    col_Z : np.array (float)
        Difference in mean densities between labeled and control gradients
    col_W_light : np.array (float)
        Mean density of all control gradients
    col_M_light : np.array (float)
        Molecular weight of DNA in control gradients
    """
    cdef int x = 1

    cdef int n = len(col_Z)
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    for i in xrange(n):
        out[i] = (col_Z[i] / col_W_light[i] + x) * col_M_light[i]
    return out


def atomFracExcess(np.ndarray[DTYPE_t, ndim=1] col_M_lab,
                   np.ndarray[DTYPE_t, ndim=1] col_M_light,
                   np.ndarray[DTYPE_t, ndim=1] col_M_heavyMax,
                   isotope='13C'):
    """Calculate atom fraction excess from qSIP data.
    
    Parameters
    ----------
    col_M_lab : np.array (float)
        Molecular weight of DNA in labeled treatment gradients
    col_M_light : np.array (float)
        Molecular weight of DNA in control treatment gradients       
    col_M_heavyMax : np.array (float)
        Theoretical molecular weight of fully-labeled DNA
    isotope : str
        13C or 18O isotope?
    """
    assert len(col_M_lab) == len(col_M_light)
    assert len(col_M_lab) == len(col_M_heavyMax)

    cdef int z = 1
    cdef double a 
    if isotope == '13C':
        a = 0.01111233
    elif isotope == '18O':
        a = 0.002000429
    else:
        msg = '"{}" Isotope not supported; only: "13C" or "18O"'
        raise TypeError, msg.format(isotope)

    cdef int n = len(col_M_lab)
    cdef double[:] out = np.empty(n, dtype=DTYPE)
    cdef unsigned int i
    cdef double x
    cdef double y
    for i in xrange(n):        
        x = col_M_lab[i] - col_M_light[i]
        y = col_M_heavyMax[i] - col_M_light[i]
        out[i] = x / y * (z - a)
    return out


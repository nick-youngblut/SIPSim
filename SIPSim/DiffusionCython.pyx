from __future__ import division
import sys
import random
import numpy as np
cimport numpy as np

#import SIPSimCpp

DTYPE = np.float
ctypedef np.float_t DTYPE_t

    
def add_diffusion_Clay(np.ndarray[DTYPE_t, ndim=2] arr, 
                       float T=298, float B=1.195e9, 
                       float G=7.87e-10, int M=882):
    """Adding diffusion to fragment buoyant density values. 
    See Clay et al., 2003 for more details.

    Paramters
    ---------
    arr : numpy array
        2d array: [[frag_buoyant_density,],[frag_length,]]
    T : float
        gradient temperature in Kelvin
    B : float
        beta coefficient
    G : float
        G coefficient
    M : float
         molecular weight per pair base pair of dry cesium DNA

    Returns
    -------
    numpy.array : [BD values that include diffusion error]
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
        diff_error = calc_diffusion_BD(arr[0,i], arr[1,i], T, B, G, M)

        # true_BD + diffusion_error_BD
        out[i] = arr[0,i] + diff_error
        
    return out


def calc_diffusion_BD(double frag_BD, double frag_len, double T, double B, double G, int M):
    """Calculating diffusion in standard deviation of buoyant_density equivalents (rho).
    Paramters
    ---------
    frag_BD : rho (buoyant density)
    frag_len : fragment length (bp)
    T : absolute temperature
    B : beta
    G : G coefficient (see Clay et al., 2003)
    M : molecular weight per base pair of dry cesium DNA

    Returns
    -------
    float : BD error due to diffusion value drawn from a normal distribution with 
            a standard deviation determined by calculated diffusion
    """    
    cdef double R = 8.3145e7
    cdef double sd_BD
    
    sd_BD = np.sqrt((frag_BD * R * T)/(np.power(B,2) * G * M * frag_len))
    return np.random.normal(0, scale=sd_BD)


def calc_sigma_Clay(double frag_BD, double frag_len,
                    double T, double B, double G, int M):
    """ Calculating diffusion error as standard deviation of buoyant_density 
    equivalents (rho). Using the Clay et al., (2003) method.
    
    Parameters
    ----------
    frag_BD = rho (buoyant density)
    frag_len = fragment length (bp)
    T = absolute temperature
    B = beta
    G = G coefficient (see Clay et al., 2003)
    M = molecular weight per base pair of dry cesium DNA
    Return:
    BD error due to diffusion value drawn from a normal distribution with a
    standard deviation determined by calculated diffusion
    """
    cdef double R = 8.3145e7
    cdef double sigma
    
    sigma = np.sqrt((frag_BD * R * T)/(B**2 * G * M * frag_len));
    return sigma


def calc_sigma_Mes(double p_p, double l, int t, double r_t, double r_b, 
                   double p_m, double B, double w):
    """ Calculating diffusion error as standard deviation of buoyant density
    equivalents (rho). Using Meselson et al., (1957) method.
    
    Parameters
    ----------
    p_p = buoyant density of the the particle at equilibrium (g ml^-1)
    l = dsDNA length (bp)
    t = time (sec)
    r_t = distance between top of gradient and center of rotation (cm)
    r_b = distance between bottom of gradient and center of rotation (cm)    
    p_m = average density of the medium (g ml^-1)
    B = beta coef. of salt forming the density gradient
    w = angular velocity 
    """
    cdef float GC
    cdef float S 
    cdef float r_c
    cdef float r_p
    cdef float sigma
    cdef float sigma_BD
    cdef float L = r_b - r_t

    GC = _BD2GC(p_p) * 100
    S = calc_S(l, GC)
    r_c = calc_R_c(r_t, r_b)
    r_p = calc_R_p(p_p, p_m, B, w, r_c)
    sigma = calc_diff_sigma(L, w, r_p, S, t, B, p_p, p_m)
    sigma_BD = sigma2BD(r_p, sigma, p_m, B, w, r_c)
    if sigma_BD < 0:
        raise ValueError, 'Sigma < 0 for BD={}, L={}'.format(p_p, l)
    return(sigma_BD)


def GC2MW(x):
    """Convert GC to molecular weight (1 bp of ds DNA)
    Parameters
    ----------
    x = Percent GC content
    """
    cdef float A = 313.2
    cdef float T = 304.2
    cdef float C = 289.2
    cdef float G = 329.2
    cdef float y = x / 100.0
    return y*(G+C) + (1-y)*(A+T)  


def GC2BD(np.ndarray[DTYPE_t, ndim=1] arr):
    """Convert G+C (% from 0-100) to buoyant density (BD)
    Equation: (GC / 100) * 0.098 + 1.66

    Parameters
    ----------
    arr = G+C values

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


def _BD2GC(double GC):
    """buoyant density to G+C (fraction). 
    Ref: Birnie and Rickwood, 1978

    Returns
    -------
    float : G+C molar fraction
    """
    cdef double x = 1.66
    cdef double y = 0.098

    return (GC - x) / y 


def calc_R_c(float r_t, float r_b):
    """Calculating the isoconcentration point.
    """
    cdef double x
    x = r_t**2 + r_t * r_b + r_b**2
    x = np.sqrt(x/3) 
    return x


def calc_R_p(float p_p, float p_m, float B, float w, float r_c):
    """Calculate the distance of the particle from the axis of rotation 
    (at equilibrium)
    """
    cdef float x
    x = ((p_p - p_m) * 2 * B / w) + r_c**2
    x = np.sqrt(x)
    return x 


def calc_S (float l, float GC):
    """Calculate the sedimentation coef based on fragment length
    """
    cdef float S
    cdef MW = GC2MW(GC)
    S = (2.8 + 0.00834 * (l * MW)**0.479) * 1e-13
    return S


def calc_diff_sigma(float L, float w, float r_p, float S, float t, 
                    float B, float p_p, float p_m):
    """
    p_p = BD
    p_m = test
    """
    cdef float nom
    cdef float denom
    cdef float x
    nom = w**2 * r_p**2 * S * t
    denom = B * abs(p_p - p_m)
    if denom == 0:
        return 0
    x = nom / denom - 1.26
    return L / np.exp(x)


def R_p2BD (float r_p, float p_m, float B, float w, float r_c):
    """Converting a distance from center of rotation of a particle to buoyant
    density (inverse of `calc_R_p`)
    """
    cdef float nom
    nom = (r_p**2 - r_c**2) * w
    return nom / (2 * B) + p_m


def sigma2BD (float r_p, float sigma, float p_m, float B, float w, float r_c):
    """Converting sigma (distance from center of rotation) to buoyant density
    """
    cdef float BD_low
    cdef float BD_high
    BD_low = R_p2BD(r_p - sigma, p_m, B, w, r_c)
    BD_high = R_p2BD(r_p + sigma, p_m, B, w, r_c)
    return BD_high - BD_low

    

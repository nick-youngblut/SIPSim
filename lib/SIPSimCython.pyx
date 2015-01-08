from __future__ import division
import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

def addIncorpBD(np.ndarray frag_BD, np.ndarray incorp_perc, double isotopeMaxBD):
    cdef int l1 = frag_BD.shape[0]
    cdef int l2 = incorp_perc.shape[0]
    cdef int i
    cdef np.ndarray res = np.zeros([l1], dtype=DTYPE)
    
    for i in xrange(l1):
        res[i] = incorp_perc[i] / 100 * isotopeMaxBD + frag_BD[i]

    return res
        

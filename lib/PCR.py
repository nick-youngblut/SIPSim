"""PCR simualation functions"""

# import
## batteries
import os,sys
import copy
import math
## 3rd party
import numpy as np
import pandas as pd
## app
import Utils


def PCR_cycle(M_n, P_n, f_0, k):
    """Increase in template concentratino after a PCR cycle.
    
    Parameters
    ----------
    M_0 : float
        The template molarity at PCR cycle 0
    P_n : float
        The primer molarity 
    f_0 : float
        The theoretical maximum PCR rxn efficiency
    k : float
        The ratio between the rate constants of reannealing and priming rxns

    Returns
    -------
    float : post-cycle template molarity 
    """
    assert M_n >= 0, 'M_n must be >= 0'
    if M_n == 0:
        return 0

    # calculating rxn efficiency
    b = k * M_n + P_n
    if b > 0:
        f_n = f_0 * P_n / b
    else:
        f_n = 0.0

    # debug
    #print (f_n, M_n, P_n, f_0, k, b)
        
    # calculating resulting template molarity
    return M_n *  math.e ** f_n
    

def run_PCR(taxa_molarities, P_0, f_0, k, n=30):
    """PCR run simulation.

    Parameters
    ----------
    taxa_molarities : list
        molarities for each taxon in community
    P_0 : float
        The inital primar molarity
    f_0 : float
        The initial PCR rxn efficiency
    k  : float 
        The ratio between the rate constants of reannealing and priming rxns
    Returns
    -------
    float : taxa molarities post-PCR
    """
    P_n = P_0 * 1e-6
    tm = list(taxa_molarities)
    for cycle in range(n):
        # debug
        #print "cycle:{}".format(cycle)

        # initial template molarity
        M_sum1 = np.sum(tm)

        # calculating new molarity for cycle (for each taxon)
        f = lambda M_0 : PCR_cycle(M_0, P_n=P_n, f_0=f_0, k=k)
        tm = [f(x) for x in tm]

        # sum of new molarity
        M_sum2 = np.sum(tm)

        # change in primer molarity (inverse of template Molarity)
        P_n = P_n - (M_sum2 - M_sum1)
        
    return tm


def PCR_sim(otu_tbl, DNA_conc_dist, DNA_conc_dist_p, primer_conc, 
            n_cycles=30, f_0=1, k=5, debug=0):
    """Simulate PCR on each gradient fraction sample.
    In-place edit of OTU table.

    Parameters
    ----------
    otu_table : OTU table object
    DNA_conc_dist : function
        numpy.random distribution used to select starting DNA 
        molarity for each sample (units = uM). 
    DNA_conc_dist_p : dict
        DNA_conc_dist parameters.
    primer_conc : float
        molarities of primers (each), units = uM.
    n_cycles : int
        number of PCR
    f_0 : float
        The theoretical maximum PCR rxn efficiency
    k : int 
        k parameter in Suzuki & Giovannoni (1996)

    Returns
    -------
    otu_table object : otu_table with edited values
    """
    # making a partial function for DNA_conc_dist
    dist_func = Utils.part_dist_func(DNA_conc_dist, DNA_conc_dist_p)

    # adding partial template molarities for each community
    f = lambda x : x * dist_func(size=1)[0] * 1e-6
    otu_tbl.apply_each_comm(f, ['rel_abund'], ['init_molarity'])

    # PCR
    f = lambda x : run_PCR(x, P_0=primer_conc, f_0=f_0, k=k, n=30)
    otu_tbl.apply_each_comm(f, ['init_molarity'], ['final_molarity'])

    # calculating new relative abundances by using final molarity as proportions
    otu_tbl.add_rel_abund(sel_index=['final_molarity'], 
                          val_index=['rel_abund'])

    # adjusting counts (absolute abundances) based on relative abundances 
    otu_tbl.adjust_abs_abund(rel_index=['rel_abund'], 
                             abs_index=['count'])

    # removing intermediate columns
    if not debug:
        cols2rm = ['init_molarity', 'final_molarity']
        otu_tbl.rm_columns(cols2rm)

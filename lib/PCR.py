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


def PCR_cycle(M_0, P_n, f_0, k):
    """Increase in template concentration after a PCR cycle.
    
    Parameters
    ----------
    M_0 : float
        The template molarity at the start of the PCR cycle
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
    assert M_0 >= 0, 'M_n must be >= 0'

    # calculating rxn efficiency
    b = k * M_0 + P_n 
    if b > 0:
        f_n = f_0 * (P_n / b) 
    else:
        f_n = 0.0

    # calculating resulting template molarity
    M_n = M_0 * math.e ** f_n 
    
    # assert
    if M_n < M_0:
        msg = 'ERROR: M_n < M_0; f_n: {}\nf_0: {}\nM_0: {}\nP_n: {}\nk: {}\nb{}' 
        sys.exit(msg.format(f_n, f_0, M_0, P_n, k, b))
    
    # return
    return M_n
    

def run_PCR(taxa_molarities, P_n, f_0, k, n_cycles=30, ratio=10):
    """PCR run simulation (PCR on 1 sample)

    Parameters
    ----------
    taxa_molarities : list
        molarities for each taxon in community
    P_n : float
        The primer molarity
    f_0 : float
        The initial PCR rxn efficiency
    k : float 
        The ratio between the rate constants of reannealing and priming rxns
    n : int
        The number of PCR cycles
    ratio : float
        Amplicon to primer length ratio.
    Returns
    -------
    float : taxa molarities post-PCR
    """
    P_n = P_n * 1e-6
    tm = list(taxa_molarities)
    for cycle in range(n_cycles):
        # initial template molarity
        M_init_sum = np.sum(tm)

        # calculating efficiency
        f_n = f_0 * (P_n / (k * M_init_sum + P_n)) *2
        
        # calculating new molarity for cycle (for each taxon)
        f = lambda M_0 : PCR_cycle(M_0, P_n=P_n, f_0=f_n, k=k)
        tm = [f(x) for x in tm]

        # sum of new molarity
        M_post_sum = np.sum(tm)

        # checking that molarity is increasing
        if M_post_sum - M_init_sum < 0:
            msg = 'Template molarity decreased through PCR cycle!'
            sys.stderr.write('ERROR: ' + msg + '\n')
            sys.exit()

        # change in primer molarity (inverse of delta template Molarity)
        P_n = P_n - (M_post_sum - M_init_sum) / ratio
        if P_n <= 0:
            msg = 'Cycle {}: Primer conc = 0 uM\n'
            sys.stderr.write(msg.format(cycle + 1))
            break

    return tm


def PCR_sim(otu_tbl, DNA_conc_dist, DNA_conc_dist_p, primer_conc, 
            n_cycles=30, f_0=1, k=5, ratio=10, debug=0):
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
    f = lambda x : run_PCR(x, P_n=primer_conc, f_0=f_0, k=k, 
                           n_cycles=n_cycles, ratio=ratio)
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

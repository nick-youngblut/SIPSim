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


def rxn_eff(f_0, P_n, M_n, k):
    """Calculating PCR rxn efficiency based on Suzuki and Giovannoni (1996).
    Args:
    f_0 -- The rxn effeciency at PCR cycle 0
    P_n -- Primer molarity at PCR cycle n
    M_n -- Template molarity at PCR cycle n
    k -- The ratio between the rate constants of reannealing and priming rxns
    Return:
    f_n -- The PCR rxn efficiency at PCR cycle n
    """ 
    try:
        f_n = f_0 * (P_n / (k * (M_n + P_n)))
    except ZeroDivisionError:
        f_n = 0.0
    return f_n


def calc_template_conc(M_0, P_0, f_0, k, n, debug=0):
    """Calculate the DNA template concentration at cycle n.
    Based on Suzuki and Giovannoni (1996).
    Args:
    M_0 -- The template molarity at PCR cycle 0
    P_n -- The primer molarity 
    f -- The PCR rxn efficiency
    k -- The ratio between the rate constants of reannealing and priming rxns
    n -- The PCR cycle number
    debug -- extra output to stderr
    Return:
    M_n -- the template molarity at cycle n
    """
    # skipping if M_0 == 0 (nothing to amplify)
    if M_0 <= 0:
        return 0

    # debug msg
    dmsg = 'cycle:{}, M_n:{}, P_n:{}, f_n:{}\n'

    # converting units (uM --> M)
    #M_0 already as uM
    P_0 = P_0 * 1e-6


    # for each PCR cycle
    M_n = M_0
    P_n = P_0
    for cycle in range(n):
        # getting efficiency
        f_n = rxn_eff(f_0, P_n, M_n, k)
        # calculating new molarity at end of cycle n
        M_n1 = M_n *  math.e ** f_n
        # calculating drop in primer conc (subtract of increase in template)
        P_n = P_n - (M_n1 - M_n)
        # setting template molarity for next cycle
        M_n = M_n1
        # debug
        if debug:
            sys.stderr.write(dmsg.format(cycle, M_n, P_n, f_n))

    return M_n



def PCR_Sim(otu_tbl, DNA_conc_dist, DNA_conc_dist_p, primer_conc, 
            n_cycles=30, f_0=1, k=5, debug=0):
    """Simulate PCR on each gradient fraction sample.
    Args:
    otu_table -- OTU table object
    DNA_conc_dist -- numpy.random distribution used to select starting DNA 
                     molarity for each sample (units = uM). 
    DNA_conc_dist_p -- DNA_conc_dist parameters.
    primer_conc -- molarities of primers (each), units = uM.
    n_cycles -- number of PCR
    f_0 -- The init PCR rxn efficiency
    k -- k parameter in Suzuki & Giovannoni (1996)
    Return:
    otu_table -- copy of otu_table with edited values
    """

    # making a partial function for DNA_conc_dist
    dist_func = Utils.part_dist_func(DNA_conc_dist, DNA_conc_dist_p)

    # adding partial template molarities for each community
    otu_tbl.add_init_molarity(dist_func)
    
    # calculating per-taxon post-PCR molarities
    f = lambda row: calc_template_conc(M_0=row['init_molarity'],
                                       P_0=primer_conc,
                                       f_0=f_0, 
                                       k=k,
                                       n=n_cycles,
                                       debug=debug)
    otu_tbl.apply_each_taxon(f, 'final_molarity')

    # calculating new relative abundances
    otu_tbl.add_rel_abund('final_molarity', val_index='rel_abund_post_PCR')
    
    # removing intermediate columns
    if not debug:
        cols2rm = ['init_molarity', 'final_molarity']
        otu_tbl.rm_columns(cols2rm)

    # writing out otu table
    otu_tbl.to_csv(sys.stdout, sep='\t', index=False)

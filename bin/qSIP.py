#!/usr/bin/env python

#--- Option parsing ---#
"""
qSIP: simulate quantitative SIP data

Usage:
  qSIP [options] <OTU_table> <OTU_subsample_table>
  qSIP -h | --help
  qSIP --version

Options:
  <OTU_table>             OTU table file ('true abundances').
  <OTU_subsample_table>   OTU table file (post-sequencing abundances;
                                           relative abundances).
  --copies=<c>            Total gene copy number.
                          [Default: 1e9]
  --reps=<r>              Number of qPCR replicates.
                          [Default: 3]
  --error_dist=<dc>       Distribution to use. 
                          (see Description for possible distributions)
                          [Default: neg_binom]
  --error_dist_p=<dp>     Distribution parameters.
                          (see Description for possible distributions)
                          [Default: alpha:0.5]
  --version               Show version.
  --debug                 Debug mode.
  -h --help               Show this screen.


Description:
  Simulate quantitative stable isotope probing (qSIP) data.

References:
  Hungate BA, Mau RL, Schwartz E, Caporaso JG, Dijkstra P, Gestel N van, et
  al. (2015). Quantitative Microbial Ecology Through Stable Isotope Probing.
  Appl Environ Microbiol AEM.02280–15.
"""

# import
## batteries
from docopt import docopt
import sys
import os
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from Utils import parseKeyValueString as distParamParse
from OTU_Table import OTU_table
import QSIP
    

def main(Uargs):
    # METHOD (step 1: apply qPCR data)
    # get total 'true' sample abundances from OTU table
    # calculate qPCR values
    ## get value(s) from neg-binom distribution
    ### `from Error_Dist import error_dist`
    #### alpha = ?
    ## [aside] write qPCR value table
    ## get mean of replicates 
    ### used to extrapolate abundances
    # transform rel-abundances 
    # write table of transformed values

    # METHOD (step 2: calculate % atom incorp)
    # NOTE: need to specify which libraries are control vs treatment
    ## comm or incorp file?
    # per taxon:
    ## Calc: total gradient copy number per taxon
    ## Calc: taxon_density = weighted average of density 
    ###      (weights=copy numbers by fraction)
    ## Calc: mean_taxon_density (summed across rep gradients)
    ## Calc: density_shift = (mean_density__treat - mean_density__control)
    ## Calc: atom_fraction_excess = see Hungate et al., AEM

    # parsing dist params
    Uargs['--error_dist_p'] = distParamParse(Uargs['--error_dist_p'])

    # error distribution function
    e_dist = error_dist(Uargs['--error_dist'],
                        Uargs['--error_dist_p'])

    # loading OTU table (s)
    otu_abs = OTU_table.from_csv(Uargs['<OTU_table>'], sep='\t')
    otu_rel = OTU_table.from_csv(Uargs['<OTU_subsample_table>'], sep='\t')

    # getting total absolute abunds for each sample
    #total_abs = otu_abs.apply_each_sample(sum, 'count')
    
    # drawing error from OTU counts
#    f = lambda x : e_dist.sample(1, x)[0]
#   total_abs_err = [f(x) for x in total_abs] 

    # Transforming relative abunds to abs abunds
    ## determining the proportional absolute abundances
    ## = total_gene_copies * taxon_rel_abundance
    
    # writing out transformed OTU (rel->abs) table
#    otu_tbl.to_csv(sys.stdout, sep='\t', index=False)


    # qSIP simulation
#    qSIP_sim()

    # writing out file
    #otu_tbl.to_csv(sys.stdout, sep='\t', index=False)
    

# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    main(Uargs)

    

        

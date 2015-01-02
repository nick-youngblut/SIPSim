#!/usr/bin/env python

#--- Option parsing ---#
"""
OTU_table: create an OTU table of gradient fractions based on simulated
fragment G+C content and isotope incorporation

Usage:
  OTU_table [options] <fragGC_file> <comm_file> <incorp_file> <frac_file>
  OTU_table -h | --help
  OTU_table --version

Options:
  <fragGC_file>       Name of file produced by fragGC subcommand.
  <comm_file>         Name of file produced by gradientComms subcommand.
  <incorp_file>       Name of file produced by isoIncorp subcommand.
  <frac_file>         Name of file produced by fractions subcommand.
  --abs_abund=<aa>    Absolute abundance of all taxa in the community. [default: 1e6]
  --g_noise=<gn>      scipy distribution function describing gradient 'noise'. [default: cauchy]
  --gn_scale=<np>     Scale parameter for the '--g_noise' distribution. [default: 0.0]
  --gc_range=<gcr>    Min-max possible G+C values post-diffusion or post-noise. [default: 0,100]
  --a_weight=<aw>     Abundance weighting for isotope incorporation.
  --isotope=<is>      Isotope incorporated by taxa (13C or 15N). [default: 13C]
  -h --help           Show this screen.
  --version           Show version.
  --debug             Debug mode

Description:

"""

# import
## batteries
from docopt import docopt
import os, sys
import logging
from collections import defaultdict
import time
## 3rd party
import numpy as np
import pandas as pd
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import SIPSim as SS

# functinos
def main(Uargs):
    
    # --abs_abund as int
    try:
        Uargs['--abs_abund'] = int(float(Uargs['--abs_abund']))
    except KeyError:
        raise KeyError('Cannot find "--abs_abund" key')
    except TypeError:
        raise TypeError('"{}" must be float-like'.format(str(Uargs['--abs_abund'])))
    
    # loading fragGC file as multivariate KDEs
    fragKDE = SS.Frag_multiKDE(Uargs['<fragGC_file>'])

    # loading community file
    comm = SS.CommTable.from_csv(Uargs['<comm_file>'], sep='\t')
    comm.set_abs_abund(Uargs['--abs_abund'])
    
    
    # loading incorp file
    incorp = SS.IsoIncorpTable.from_csv(Uargs['<incorp_file>'], sep='\t')    
    
    # loading fraction file
    frac = SS.FracTable.from_csv(Uargs['<frac_file>'], sep='\t')
    
    # initializing OTU table class
    OTU = SS.OTU_table(frac,
                    g_noise=Uargs['--g_noise'],
                    gn_scale=Uargs['--gn_scale'],
                    gn_range=Uargs['--gn_range'],
                    abund_weight=Uargs['--a_weight'],                    
                    isotope=Uargs['--isotope'])
    
    # checking on library overlap
    if not OTU.checkLibOverlap([
            [x for x in comm.iter_libraries()],
            [x for x in incorp.iter_libraries()],
            [x for x in frac.iter_libraries()]]):
        logging.warning('Not all tables contain the same library IDs')

        
    # iter by library:
    isotopeMaxBD = OTU.get_isotopeMaxBD()
    OTU_counts = dict()
    for libID in comm.iter_libraries():
        sys.stderr.write('Processing library: "{}"\n'.format(libID))
        
        # make dataframe for OTU counts: taxa X fractions
        OTU_counts[libID] = OTU.make_emptyCountTable(comm.get_unique_taxon_names(),
                                                     frac.get_libFracIDs(libID))

        # all values to integers
        OTU_counts[libID] = OTU_counts[libID].astype(int)

        # iter by taxon:
        for taxon_name in comm.iter_taxa(libID):
            sys.stderr.write('  Processing taxon: "{}"\n'.format(taxon_name))
            
            taxonAbsAbund = comm.get_taxonAbund(libID, taxon_name)

            # sampling fragment GC & length values from taxon-specific KDE
            t0 = time.time()
            GC_len_arr = fragKDE.sampleTaxonKDE(taxon_name, size=taxonAbsAbund)            
            
            # sampling intra-taxon incorp for taxon; return: iterator
            t1 = time.time()
            incorp_vals = incorp.sample_incorpFunc(libID, taxon_name, n_samples=taxonAbsAbund)
            
            # iter GC value:
            t2 = time.time()
            for (frag_gc,frag_len) in GC_len_arr:
                # simulating diffusion on GC
                frag_gc = OTU.add_diffusion(loc=0, frag_len=frag_len)
#                frag_gc += np.random.normal(loc=0, scale=44500/frag_len)                

                # simulating noise
                frag_gc = OTU.sample_g_noise_func(frag_gc, loc=0)

                print frag_gc; sys.exit()
                
                # raw BD based on GC
                BD = frag_gc / 100 * 0.098 + 1.66
                
                # BD + BD shift from isotope incorporation
                ## TODO: implement abundance-weighting
                incorp_perc = incorp_vals.next()
                BD = BD + isotopeMaxBD * (incorp_perc / 100)

                
                # simulate diffusion
                #BD += np.random.normal(loc=0, scale=44500/frag_len) / 100 * 0.098 + 1.66

#                print BD
                
                # simulate noise
                BD += OTU.sample_g_noise_func(0)[0]

 #               print BD;
                sys.exit()
                
                # simulating gradient noise
                #BD = OTU.sample_g_noise_func(BD)[0]
                
                # determine the fraction that would contain the fragment
                #fracID = frac.which_frac(libID, BD)

                # adding to OTU count table
                try:
                    OTU_counts[libID][fracID].loc[taxon_name] += 1
                except ValueError:
                    pass

            t3 = time.time()
            #print [t1 - t0, t2 - t1, t3 - t2]; sys.exit()
                    
        # replace fraction IDs in column names with BD_min-BD_max
        OTU_counts[libID].columns = frac.fracID2BDminmax(libID, OTU_counts[libID].columns)

                    
    # combine the dataframes and write
    pd.concat(OTU_counts.values(), axis=1).to_csv(sys.stdout, sep='\t')


        
# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    main(Uargs)
    

        
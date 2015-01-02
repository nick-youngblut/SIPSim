# import
## batteries
import os,sys
import math
import logging
import time
## 3rd party
import pandas as pd
import numpy as np
## application
import FragGC
from CommTable import CommTable
from FracTable import FracTable
from IsoIncorpTable import IsoIncorpTable
from OTU import OTU_table


# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


# main
def main(Uargs):    
    # --abs_abund as int
    try:
        Uargs['--abs_abund'] = int(float(Uargs['--abs_abund']))
    except KeyError:
        raise KeyError('Cannot find "--abs_abund" key')
    except TypeError:
        raise TypeError('"{}" must be float-like'.format(str(Uargs['--abs_abund'])))
    
    # loading fragGC file as multivariate KDEs
    fragKDE = FragGC.Frag_multiKDE(Uargs['<fragGC_file>'])

    # loading community file
    comm = CommTable.from_csv(Uargs['<comm_file>'], sep='\t')
    comm.set_abs_abund(Uargs['--abs_abund'])
    
    
    # loading incorp file
    incorp = IsoIncorpTable.from_csv(Uargs['<incorp_file>'], sep='\t')    
    
    # loading fraction file
    frac = FracTable.from_csv(Uargs['<frac_file>'], sep='\t')
    
    # initializing OTU table class
    OTU = OTU_table(frac,
                    g_noise=Uargs['--g_noise'],
                    gn_scale=Uargs['--gn_scale'],
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
                frag_gc = OTU.add_diffusion(frag_gc, frag_len, loc=0)

                # simulating noise (if needed)
                if OTU.gn_scale_nonzero:
                    frag_gc = OTU.sample_g_noise_func(frag_gc, loc=0)
                
                # raw BD based on GC
                BD = frag_gc / 100 * 0.098 + 1.66
                
                # BD + BD shift from isotope incorporation
                ## TODO: implement abundance-weighting
                incorp_perc = incorp_vals.next()
                BD = BD + isotopeMaxBD * (incorp_perc / 100)
                            
                # determine the fraction that would contain the fragment
                fracID = frac.which_frac(libID, BD)

                # adding to OTU count table
                try:
                    OTU_counts[libID][fracID].loc[taxon_name] += 1
                except ValueError:
                    logging.warning('BD value {} does not fall into any fraction'.format(BD))
                    pass

            t3 = time.time()
            #print [t1 - t0, t2 - t1, t3 - t2]; sys.exit()
                    
        # replace fraction IDs in column names with BD_min-BD_max
        OTU_counts[libID].columns = frac.fracID2BDminmax(libID, OTU_counts[libID].columns)

                    
    # combine the dataframes and write
    pd.concat(OTU_counts.values(), axis=1).to_csv(sys.stdout, sep='\t')



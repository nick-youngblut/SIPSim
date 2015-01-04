# import
## batteries
import os,sys
import logging
import time
from collections import defaultdict, Counter
from itertools import starmap, imap
## 3rd party
import pandas as pd
import numpy as np
import parmap
from numba import jit, vectorize
## application
import FragGC
from CommTable import CommTable
from FracTable import FracTable
from IsoIncorpTable import IsoIncorpTable
import OTU


# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

# supporting functions
def makeEmptyCountTable(comm, frac):
    """Making a dict of 2d np arrays to hold the taxon count data.
    All counts = 0. Dimensions = library : [taxon][fraction]
    Args:
    comm -- CommTable instance
    frac -- FracTable instance
    Return:
    2d-array of zeros (n-taxon, n-fraction)
    """
    # assertions
    assert isinstance(comm, CommTable), 'ERROR comm is not a CommTable instance'
    assert isinstance(frac, FracTable), 'ERROR comm is not a FracTable instance'
    
    # making count table
    taxon_name_len = len(comm.get_unique_taxon_names())
    countTable = dict()
    for libID in comm.iter_libraries():
        fracID_len = len(frac.get_libFracIDs(libID))
        countTable[libID] = np.zeros((taxon_name_len, fracID_len))
        
    return countTable

    
def calcFragBD_OLD():
    # simulating diffusion on GC
    frag_gc = OTU.add_diffusion(frag_gc, frag_len, loc=0)
    
    # simulating noise (if needed)
    if OTU.gn_scale_nonzero:
        frag_gc = OTU.sample_g_noise_func(frag_gc, loc=0)
        
    # raw BD based on GC
    BD = frag_gc / 100.0 * 0.098 + 1.66
        
    # BD + BD shift from isotope incorporation
    ## TODO: implement abundance-weighting
    incorp_perc = incorp_vals.next()
    BD = BD + isotopeMaxBD * (incorp_perc / 100)
    
    # determine the fraction that would contain the fragment
    fracID = frac.which_frac(libID, BD)
                
    # adding to OTU count table
    try:
        OTU_counts[libID][taxon_name_idx,fracID-1] += 1
    except ValueError:
        logging.warning('BD value {} does not fall into any fraction'.format(BD))
        pass

        
@jit('f8[:](f8[:])')
def GC2BD(frag_gc):
    """Conversion of G+C (0-100) to buoyant density.
    Args:
    frag_gc -- numpy array of fragment G+C values
    """
    # TODO: jit not speeding up calculation, need to fix
    return frag_gc / 100.0 * 0.098 + 1.66

    
@jit('f8[:](f8[:],f8[:],f8)')
def addIncorpBD(frag_BD, incorp_perc, isotopeMaxBD):
    """Adding BD from isotope incorporation to 0%-incorp BD values.
    Args:
    frag_BD -- array of frag_BD values
    incorp_perc -- array of isotope incorp percentages
    isotopeMaxBD -- float of max BD possible for isotope
    Return:
    1d array of updated BD values
    """
    n = frag_BD.shape[0]
    result = np.zeros([n])
    for i in xrange(n):
        result[i] = frag_BD[i] + isotopeMaxBD * (incorp_perc[i] / 100)
    
    return result


def binNum2ID(frag_BD_bins, libFracBins):
    """Convert Counter(np.digitize()) for frag_BD  to the fraction BD-min-max.
    Args:
    frag_BD_bins -- dict of counts from np.digitize() binning
    libFracBins -- ordered set of fraction start-ends
    Return:
    dict of {fractionID : fragment_counts}
    """
    return {'{}-{}'.format(libFracBins[k-1],libFracBins[k]):v for (k,v) in frag_BD_bins.items()}
    
    

# main
#@profile
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
    OTUt = OTU.OTU_table(frac,
                         g_noise=Uargs['--g_noise'],
                         gn_scale=Uargs['--gn_scale'],
                         abund_weight=Uargs['--a_weight'],                    
                         isotope=Uargs['--isotope'])
    
    # checking on library overlap
    if not OTUt.checkLibOverlap([
            [x for x in comm.iter_libraries()],
            [x for x in incorp.iter_libraries()],
            [x for x in frac.iter_libraries()]]):
        logging.warning('Not all tables contain the same library IDs')
    
        
    # iter by library:
    isotopeMaxBD = OTUt.get_isotopeMaxBD()
    u_taxon_names = comm.get_unique_taxon_names()
    OTU_counts = makeEmptyCountTable(comm, frac)
    
    for libID in comm.iter_libraries():
        sys.stderr.write('Processing library: "{}"\n'.format(libID))

        # fraction bin list for library
        libFracBins = [x for x in frac.BD_bins(libID)]
        
        # all values to integers
        #OTU_counts[libID] = OTU_counts[libID].astype(int)

        # iter by taxon:
        for (taxon_name_idx,taxon_name) in enumerate(u_taxon_names):
            sys.stderr.write('  Processing taxon: "{}"\n'.format(taxon_name))
            
            taxonAbsAbund = comm.get_taxonAbund(libID, taxon_name)

            # sampling fragment GC & length values from taxon-specific KDE
            t0 = time.time()
            GC_len_arr = fragKDE.sampleTaxonKDE(taxon_name, size=taxonAbsAbund)            
            
            # sampling intra-taxon incorp for taxon; return: iterator
            t1 = time.time()
            incorp_val_iter = incorp.sample_incorpFunc(libID, taxon_name, n_samples=taxonAbsAbund)

            # GC --> BD
            t2 = time.time()

            # simulating diffusion
#            frag_gc = np.concatenate([x for x in starmap(OTU.add_diffusion,
#                                                       GC_len_arr.tolist())])            
            frag_gc = np.concatenate(parmap.starmap(OTU.add_diffusion,
                                                    GC_len_arr.tolist(),
                                                    processes=int(Uargs['--threads']),
                                                    chunksize=10000))

            t3 = time.time()
            
            # simulating general noise in the gradient column
            if OTUt.gn_scale_nonzero:
                frag_gc = map(OTU.sample_g_noise_func, frag_gc)
                
            # GC to BD values
            frag_BD = GC2BD(frag_gc)

            # BD + BD shift from isotope incorporation
            ## TODO: implement abundance-weighting
            t4 = time.time()            

            incorp_perc = np.array(list(incorp_val_iter))
            frag_BD = addIncorpBD(frag_BD, incorp_perc, isotopeMaxBD)

            t5 = time.time()
            
            # group by fraction
            frag_BD_bins = Counter(np.digitize(frag_BD, libFracBins))
            frag_BD_bins = binNum2ID(frag_BD_bins, libFracBins)            

            
            t6 = time.time()
            print [t1 - t0, t2 - t1, t3 - t2, t4 - t3, t5 - t4, t6 - t5]; #sys.exit()
#            print frag_BD; sys.exit()
                
            # TODO:
            ## convert to pd.DataFrame; combine; write out

            
                                
        # replace fraction IDs in column names with BD_min-BD_max
#        OTU_counts[libID] = pd.DataFrame(OTU_counts[libID])
#        OTU_counts[libID].index = u_taxon_names
#        OTU_counts[libID].columns = frac.fracID2BDminmax(libID)

                    
    # combine the dataframes and write
#    pd.concat(OTU_counts.values(), axis=1).to_csv(sys.stdout, sep='\t')




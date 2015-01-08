# import
## batteries
import os,sys
import logging
import time
from glob import glob
from collections import defaultdict, Counter
from itertools import imap, izip
from functools import partial
## 3rd party
import pandas as pd
import numpy as np
import multiprocessing as mp
## application
import FragGC
from CommTable import CommTable
from FracTable import FracTable
from IsoIncorpTable import IsoIncorpTable
import OTU

import RandCpp


# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


# supporting functions
def makeEmptyCountTable(comm, frac):  # depreciated
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

        
def GC2BD(frag_gc):
    """Conversion of G+C (0-100) to buoyant density.
    Args:
    frag_gc -- numpy array of fragment G+C values
    """
    # TODO: jit not speeding up calculation, need to fix
    return frag_gc / 100.0 * 0.098 + 1.66

    
def addIncorpBD(frag_BD_incorp, isotopeMaxBD):
    """Adding BD from isotope incorporation to 0%-incorp BD values.
    Args:
    frag_BD -- array of frag_BD values
    incorp_perc -- array of isotope incorp percentages
    isotopeMaxBD -- float of max BD possible for isotope
    Return:
    1d array of updated BD values
    """
    #n = frag_BD.shape[0]
    #result = np.zeros([n])
    #for i in xrange(n):
    #    result[i] = frag_BD[i] + isotopeMaxBD * (incorp_perc[i] / 100)
    return frag_BD_incorp[0] + isotopeMaxBD * (frag_BD_incorp[1] / 100)


def binNum2ID(frag_BD_bins, libFracBins):
    """Convert Counter(np.digitize()) for frag_BD  to the fraction BD-min-max.
    Args:
    frag_BD_bins -- dict of counts from np.digitize() binning
    libFracBins -- ordered set of fraction start-ends
    Return:
    dict of {fractionID : fragment_counts}
    """
    return {'{0:.3f}-{1:.3f}'.format(libFracBins[k-1],libFracBins[k]):v for (k,v) in frag_BD_bins.items()}
    
    

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


    # parallel
    chunkSize = int(Uargs['--chunksize'])
    mpPool= mp.Pool(processes=int(Uargs['--threads']))
        
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
    OTU_counts = []  # list of all library-specific OTU_count dataframes
    
    for libID in comm.iter_libraries():
        sys.stderr.write('Processing library: "{}"\n'.format(libID))

        # fraction bin list for library
        libFracBins = [x for x in frac.BD_bins(libID)]
        
        # creating a dataframe of fractions
        fracBins = ['{0:.3f}-{1:.3f}'.format(libFracBins[i-1],libFracBins[i]) for i in xrange(len(libFracBins))]
        
        lib_OTU_counts = pd.DataFrame({'fractions':fracBins})
         
        
        # iter by taxon:
        for (taxon_name_idx,taxon_name) in enumerate(u_taxon_names):
            t_start = time.time()
            sys.stderr.write('  Processing taxon: "{}"\n'.format(taxon_name))
                        
            taxonAbsAbund = comm.get_taxonAbund(libID, taxon_name)
            sys.stderr.write('    N-fragments:   {}\n'.format(taxonAbsAbund))
            
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
            #frag_gc = np.concatenate(parmap.map(OTU.add_diffusion,
            #                                    GC_len_arr,
            #                                    processes=int(Uargs['--threads']),
            #                                    chunksize=10000))
            #frag_gc = parmap.map(OTU.add_diffusion,
            #                     GC_len_arr,
            #                     parallel=True,
            #                     processes=int(Uargs['--threads']),
            #                     chunksize=100)
            frag_gc = mpPool.imap(OTU.add_diffusion,
                                  GC_len_arr,
                                  chunksize=chunkSize)
            
            t3 = time.time()
            
            # simulating general noise in the gradient column
            if OTUt.gn_scale_nonzero:
                frag_gc = mpPool.imap(OTU.sample_g_noise_func,
                                      frag_gc,
                                      chunksize=chunkSize)

                
            # GC to BD values
           # frag_BD = GC2BD(frag_gc)
            frag_BD = mpPool.imap(GC2BD, 
                                  frag_gc,
                                  chunksize=chunkSize)

            t4 = time.time()            
            
            # BD + BD shift from isotope incorporation
            ## TODO: implement abundance-weighting
            #incorp_perc = np.array(list(incorp_val_iter))
#            frag_BD = addIncorpBD(frag_BD, incorp_perc, isotopeMaxBD)
            addIncorpBD_p = partial(addIncorpBD,  isotopeMaxBD=isotopeMaxBD)

            frag_BD = mpPool.imap(addIncorpBD_p,
                                  izip(frag_BD, incorp_val_iter),
                                  chunksize=chunkSize)

            t5 = time.time()
            
            # group by fraction
            frag_BD = np.concatenate([x for x in frag_BD])
            frag_BD_bins = Counter(np.digitize(frag_BD, libFracBins))
            frag_BD_bins = binNum2ID(frag_BD_bins, libFracBins)            
            
            t6 = time.time()
            print [t1 - t0, t2 - t1, t3 - t2, t4 - t3, t5 - t4, t6 - t5]; #sys.exit()

            # converting to a pandas dataframe
            frag_BD_bins = frag_BD_bins.items()
            df = pd.DataFrame(frag_BD_bins)
            df.columns = ['fractions',taxon_name]
            df.iloc[:,1] = df.applymap(str).iloc[:,1]   # must convert values to dtype=object
            
            lib_OTU_counts = pd.merge(lib_OTU_counts, df, how='outer', on='fractions') 

            # status
            t_end = time.time()
            sys.stderr.write('    Time elapsed:  {0:.1f} sec\n'.format(t_end - t_start))
            
            

        # formatting completed dataframe of OTU counts for the library
        lib_OTU_counts.fillna(0, inplace=True)
        lib_OTU_counts = pd.melt(lib_OTU_counts, id_vars=['fractions'])
        lib_OTU_counts.columns = ['fractions', 'taxon', 'count']
        lib_OTU_counts['library'] = libID
        lib_OTU_counts = lib_OTU_counts[['library','fractions','taxon','count']]
        lib_OTU_counts.sort(['fractions'])
        ## removing end_fract-start_frac bin
        lib_OTU_counts = lib_OTU_counts.iloc[1:,:]

        # making dict of library-specific dataframes
        OTU_counts.append(lib_OTU_counts)

        
    # combining library-specific dataframes and writing out long form of table
    ## TODO: option for wide form of table
    pd.concat(OTU_counts, ignore_index=False).to_csv(sys.stdout, sep='\t', index=False)
            




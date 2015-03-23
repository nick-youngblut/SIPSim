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
import numexpr as ne
import multiprocessing as mp
from intervaltree import Interval, IntervalTree
## application
import FragGC
from CommTable import CommTable
from FracTable import FracTable
from IsoIncorpTable import IsoIncorpTable
import OTU
import SIPSimCpp
#import SIPSimCython

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


def make_diffusion_dists(start=1,end=60001,step=100,maxRange=2000000):
    """Making interval tree of normal distributions (loc=0, sigma=*see below).
    Will be used to model the observed G+C variance caused by diffusion.
    The interval tree is to limit the number of pre-generated distributions that
    will need to be created.
    Using equation: sigma = 44500 / L, where L = fragment length.
    Args:
    start -- start of sequence
    end -- end of sequence
    step -- sequence step
    maxRange -- last (largest) range added to interval tree for any outlier large fragments
    Return:
    interval tree object
    
    """
    # initialize itree
    itree = IntervalTree()
    
    for i in xrange(start,end,step):
        itree[i:i+step] = partial(np.random.normal, loc=0, scale=44500.0/i)

    # max itree
    itree[end+1:maxRange] = partial(np.random.normal, loc=0, scale=44500.0/end)

    return itree
    
    
    
    
def binNum2ID(frag_BD_bins, libFracBins):
    """Convert Counter(np.digitize()) for frag_BD  to the fraction BD-min-max.
    Args:
    frag_BD_bins -- dict of counts from np.digitize() binning
    libFracBins -- ordered set of fraction start-ends
    Return:
    dict of {fractionID : fragment_counts}
    """
    msg = '{0:.3f}-{1:.3f}'
    return {msg.format(libFracBins[k-1],libFracBins[k]):v for (k,v) in frag_BD_bins.items()}



    
#--- main ---#
@profile
def main(Uargs):
    # --abs_abund as int
    try:
        Uargs['--abs_abund'] = int(float(Uargs['--abs_abund']))
    except KeyError:
        raise KeyError('Cannot find "--abs_abund" key')
    except TypeError:
        raise TypeError('"{}" must be float-like'.format(Uargs['--abs_abund']))


    # fragment info log file
    if Uargs['--log'] != 'None':
        logfh = open(Uargs['--log'], 'w')
        logfh.write("\t".join(['library','taxon','fragment_GC','fragment_length']))
        logfh.write("\n")
    else:
        logfh = None
        
    # parallel
    chunkSize = int(Uargs['--chunksize'])
    mpPool= mp.Pool(processes=int(Uargs['--threads']))
    
    # loading community file
    sys.stderr.write('Loading files...\n')
    comm = CommTable.from_csv(Uargs['<comm_file>'], sep='\t')
    comm.set_abs_abund(Uargs['--abs_abund'])
        
    # loading incorp file
    incorp = IsoIncorpTable.from_csv(Uargs['<incorp_file>'], sep='\t')    
    
    # loading fraction file
    frac = FracTable.from_csv(Uargs['<frac_file>'], sep='\t')
    
    # loading fragGC file as multivariate KDEs
    sys.stderr.write('Creating 2d-KDEs of fragment GC & length...\n')
    fragKDE = FragGC.Frag_multiKDE(Uargs['<fragGC_file>'])
    
    # creating interval tree of normal distributions to model diffusion
    diffusionDists = make_diffusion_dists()
    
    # initializing OTU table class
    OTUsim = OTU.OTU_sim(frac,
                         g_noise=Uargs['--g_noise'],
                         gn_scale=Uargs['--gn_scale'],
                         abund_weight=Uargs['--a_weight'],                    
                         isotope=Uargs['--isotope'])
    
    # checking on library overlap
    if not OTUsim.checkLibOverlap([
            [x for x in comm.iter_libraries()],
            [x for x in incorp.iter_libraries()],
            [x for x in frac.iter_libraries()]]):
        logging.warning('Not all tables contain the same library IDs')
    
        
    # iter by library:
    sys.stderr.write('Creating OTUs...\n')
    isotopeMaxBD = OTUsim.get_isotopeMaxBD()
    u_taxon_names = comm.get_unique_taxon_names()
    OTU_counts = []  # list of all library-specific OTU_count dataframes
    
    for libID in comm.iter_libraries():
        sys.stderr.write('Processing library: "{}"\n'.format(libID))

        # fraction bin list for library
        libFracBins = [x for x in frac.BD_bins(libID)]
        
        # creating a dataframe of fraction bins
        func = lambda x: '{0:.3f}-{1:.3f}'.format(libFracBins[x-1],libFracBins[x])
        fracBins = [func(i) for i in xrange(len(libFracBins))][1:]        
        lib_OTU_counts = pd.DataFrame({'fractions':fracBins})
        
        # iter by taxon:
        for (taxon_name_idx,taxon_name) in enumerate(u_taxon_names):
            t_start = time.time()
            sys.stderr.write('  Processing taxon: "{}"\n'.format(taxon_name))
                        
            taxonAbsAbund = comm.get_taxonAbund(libID, taxon_name)
            sys.stderr.write('    N-fragments:   {}\n'.format(taxonAbsAbund))
                
            # sampling fragment GC & length values from taxon-specific KDE
            try:
                GC_len_arr = fragKDE.sampleTaxonKDE(taxon_name, size=taxonAbsAbund)
            except ValueError:
                GC_len_arr = None
                
            # calc BD 
            if GC_len_arr is None or GC_len_arr.shape[1] == 0:
                frag_BD = np.zeros(1)
            else:
                ## logging if needed
                if logfh is not None:
                    logfh.write("\n".join( ["\t".join([libID, taxon_name, str(x),str(y)]) \
                                            for x,y in np.transpose(GC_len_arr)] ))
                    logfh.write("\n")
                
                # simulating diffusion; calc BD from frag GC
                f = lambda x: SIPSimCpp.add_diffusion(x[0], x[1])
                frag_BD = np.apply_along_axis(f, 0, GC_len_arr) / 100.0 * 0.098 + 1.66
                GC_len_arr = ()
            
                # BD + BD shift from isotope incorporation
                ## TODO: implement abundance-weighting            
                incorp_vals = np.ravel(incorp.sample_incorpFunc(libID,
                                                       taxon_name,
                                                       n_samples=taxonAbsAbund))
                        
                frag_BD += incorp_vals / 100.0 * isotopeMaxBD

                
            # group by fraction
            frag_BD_bins = Counter(np.digitize(frag_BD, libFracBins))
            frag_BD = ()
            frag_BD_bins = binNum2ID(frag_BD_bins, libFracBins)            
            
            # converting to a pandas dataframe
            frag_BD_bins = frag_BD_bins.items()
            df = pd.DataFrame(frag_BD_bins)
            df.columns = ['fractions',taxon_name]
            
            # adding to dataframe
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

        # making dict of library-specific dataframes
        OTU_counts.append(lib_OTU_counts)

        
    # combining library-specific dataframes and writing out long form of table
#    pd.concat(OTU_counts, ignore_index=False).to_csv(sys.stdout, sep='\t', index=False)
            

    # finish up
    if logfh is not None:        
        logfh.close()
        sys.stderr.write( 'Fragment info written to: {}\n'.format(logfh.name))



# import
## batteries
import os,sys
import logging
import time
from glob import glob
from collections import defaultdict, Counter
from functools import partial
## 3rd party
import pandas as pd
import numpy as np
#import numexpr as ne
#import multiprocessing as mp
#from intervaltree import Interval, IntervalTree
## application
import FragGC
from CommTable import CommTable
from FracTable import FracTable
from IsoIncorpTable import IsoIncorpTable
import OTU
import SIPSimCpp
import SIPSimCython

# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

    
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


def status(msg, startTime):
    msgs = {'kde':'GC/fragment_length KDE sampled',
            'diffusion':'diffusion added to BD values',
            'incorp':'isotope incorporation added to BD values',
            'bin':'binned BD values',
            'final':'taxon finished'}

    nowTime = time.time()
    timeDiff = '{0:.1f}'.format(nowTime - startTime)
    
    try:
        x = '     Elapsed: {0:>7} sec => {1}\n'
        sys.stderr.write(x.format(timeDiff, msgs[msg.lower()]))
    except KeyError:
        raise KeyError('Cannot find status for key "{}"'.format(msg))

    return nowTime
        


    
#--- main ---#
#@profile
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
            status('KDE', t_start)
                
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
#                f = lambda x: SIPSimCpp.add_diffusion(x[0], x[1])
#               frag_BD = np.apply_along_axis(f, 0, GC_len_arr) / 100.0 * 0.098 + 1.66
                    
                frag_BD = SIPSimCython.add_diffusion_wrapper(GC_len_arr)
                GC_len_arr = ()
                status('diffusion', t_start)
                
                # BD + BD shift from isotope incorporation
                ## TODO: implement abundance-weighting
                incorp_vals = np.ravel(incorp.sample_incorpFunc(libID,
                                                       taxon_name,
                                                       n_samples=taxonAbsAbund))
                        
                frag_BD += incorp_vals / 100.0 * isotopeMaxBD
                status('incorp', t_start)
                incorp_vals = ()
                
                
            # group by fraction
            frag_BD_bins = Counter(np.digitize(frag_BD, libFracBins))
            frag_BD = ()
            frag_BD_bins = binNum2ID(frag_BD_bins, libFracBins)
            status('bin', t_start)
            
            # converting to a pandas dataframe
            frag_BD_bins = frag_BD_bins.items()
            df = pd.DataFrame(frag_BD_bins)
            df.columns = ['fractions',taxon_name]
            
            # adding to dataframe
            df.iloc[:,1] = df.applymap(str).iloc[:,1]   # must convert values to dtype=object
            lib_OTU_counts = pd.merge(lib_OTU_counts, df, how='outer', on='fractions')

            # end of loopo
            status('final', t_start)            
            

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



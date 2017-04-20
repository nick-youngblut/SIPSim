#!/usr/bin/env python

# import
## batteries
import sys,os
import cPickle as pickle 
from functools import partial
## 3rd party
import scipy.stats as stats
import numpy as np
import dill
from pathos.multiprocessing import ProcessingPool
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import DiffusionCython as DC
import Utils

import numpy as np
from itertools import cycle, izip
import collections

# functions
    
def make_kde(x, diff_index, BD_bins, len_bins, n, bw_method):
    """Make a scipy.stats.gaussian_kde from BD values.
    
    Parameters
    ----------
    x : list
        [taxon_name, kde]
        taxon_name -- taxon name
        kde -- scipy.stats kde object of BD values or None
    diff_index : dict{dict}
        Index of gaussian error distributions for all BD/length bins
    n : int
        sample size
    bw_method : str or function
         KDE bandwidth method (see scipy.stats)

    Returns
    -------
    k : tuple 
        (taxon_name, kde object of BD values)
    """
    # input unpacking a type checking
    try:
        taxon_name,kde = x
    except ValueError:
        msg = '"x" must be (taxon_name, kde)'
        raise ValueError, msg
    try:
        bw_method = float(bw_method)
    except (TypeError, ValueError) as e:
        pass 
    try: 
        bw_method = bw_method.lower()
    except AttributeError:
        pass

    # loading in diffusion error index
    with open(diff_index, 'rb') as inFH:
        diff_index = dill.load(inFH)

    # status
    msg = 'Processing: {}\n'
    sys.stderr.write(msg.format(taxon_name))
    if kde is None:
        return (taxon_name, None)
     
    # making new kde
    kdeBD = stats.gaussian_kde(add_diffusion(kde.resample(size=n), 
                                             diff_index, 
                                             BD_bins, len_bins),
                               bw_method=bw_method)

    return (taxon_name, kdeBD)


def add_diffusion(arr, diff_index, BD_bins, len_bins):    
    BD_val_bins =  np.digitize(arr[0], bins=BD_bins)
    len_val_bins = np.digitize(arr[1], bins=len_bins)

    # applying error
    for i in xrange(len(BD_val_bins)):
        BD_i = BD_val_bins[i]
        len_i = len_val_bins[i]
        try:
            # adding BD error drawn from normal distribution to BD value
            err = diff_index[BD_i][len_i](size=1)[0]
            arr[0][i] = arr[0][i] + err
        except KeyError:
            raise KeyError, 'Cannot find {}:{}'.format(BD_i, len_i)
    return arr[0]


def create_diff_index(BD_bins, len_bins, method, 
                      B, D, w, r_min, r_max, t, T, G, M):
    """Creating an index of normal distributions for determining error.
    Apply error distribution to all fragments that fall into that particular
    BD+length bin. Each normal has a mean of 0; just sigma varies.
    
    Parameters
    ----------
    BD_bins : iterable
        Bin start values for fragment BD bins
    len_bins : iterable
        Bin start values for fragment length bins
    method : str
        Method to calculate diffusion error. Either Meselson or Clay.
    
    """
    diff_index = collections.defaultdict(dict)
    index_size = 0
    for BD_i in xrange(len(BD_bins)):
        for L_i in xrange(len(len_bins)):
            BD = BD_bins[BD_i]
            L = len_bins[L_i]
            
            # detemine sigma based on on BD & ultra-cfg conditions
            if method.lower().startswith('clay'):
                sigma = DC.calc_sigma_Clay(BD, L, T, B, G, M)
            elif method.lower().startswith('mes'):
                sigma = DC.calc_sigma_Mes(p_p=BD, l=L, t=t, r_t=r_min, 
                                          r_b=r_max, p_m=D, B=B, w=w)
            else:
                raise ValueError, 'Do not know method: {}'.format(method)

            if sigma == 0:
                sigma = 1e-10
            diff_index[BD_i][L_i] = partial(np.random.normal, loc=0, 
                                            scale=sigma)
            index_size += 1
    # status
    sys.stderr.write('Index size: {}\n'.format(index_size))
        
    return diff_index


def parse_bin_range(x):
    x = x.split(',')
    if len(x) != 3:
        raise ValueError, 'The format must be "start,stop,step"'
    x = [float(y) for y in x]
    return x


def main(args):    
    """Main function for adding diffusion error to 
    a KDE of buoyant density values.

    Parameters:
    args : dict
        See ``diffusion`` subcommand
    """
    kde2d = Utils.load_kde(args['<fragment_kde>'])

    # creating a diffusion index of guassian distributions    
    start,stop,step = parse_bin_range(args['--BD_range'])
    BD_bins = np.arange(start, stop, step)
    start,stop,step = parse_bin_range(args['--len_range'])
    len_bins = np.arange(start, stop, step)
    diff_index = create_diff_index(BD_bins, len_bins,  
                                   method=args['-m'], 
                                   B=float(args['-B']), 
                                   D=float(args['-D']),
                                   w=float(args['-w']),
                                   r_min=float(args['--r_min']),
                                   r_max=float(args['--r_max']),
                                   t=int(args['-t']),
                                   T=float(args['-T']), 
                                   G=float(args['-G']), 
                                   M=float(args['-M']))

    diff_index_file = 'diffusion_index.txt'
    with open(diff_index_file, 'wb') as outFH:
        dill.dump(diff_index, outFH)


    # difussion calc in parallel
    pfunc = partial(make_kde, 
                    diff_index=diff_index_file,
                    BD_bins=BD_bins,
                    len_bins=len_bins,
                    n = int(args['-n']),
                    bw_method=args['--bw'])
    ## pool
    pool = ProcessingPool(nodes=int(args['--np']))
    if args['--debug']:
        KDE_BD = map(pfunc, kde2d.items())
    else:
        KDE_BD = pool.map(pfunc, kde2d.items())

    # pickling output
    dill.dump({taxon:KDE for taxon,KDE in KDE_BD}, sys.stdout)
    


        
if __name__ == '__main__':
#    T -- gradient temperature in Kelvin
#    B -- beta coefficient
#    G -- G coefficient
#    M -- molecular weight per pair base pair of dry cesium DNA

    frag_GC = np.array([50.0] * 10)
    frag_BD = SSC.GC2BD(frag_GC)

    def add_frag_len(frag_BD, frag_len):
        frag_lens = np.array([frag_len] * len(frag_BD))
        return np.vstack((frag_BD,frag_lens))
                          
    fragL1 = add_frag_len(frag_BD, 1000)
    fragL4 = add_frag_len(frag_BD, 4000)
    fragL44 = add_frag_len(frag_BD, 44500)
    fragL100 = add_frag_len(frag_BD, 100000)
 
    print "--BD standard deviation--"
    msg = 'Frag-len: {}, BD s.d.: {}'
    print msg.format(1000, np.std(fragL1[0]))
    print msg.format(4000, np.std(fragL4[0]))
    print msg.format(44500, np.std(fragL44[0]))
    print msg.format(100000, np.std(fragL100[0]))

    T = 298
    B = 1.195e9
    G = 7.87e-10
    M = 882
                        
    fragL1_dif = SSC.add_diffusion(fragL1, T, B, G, M)
    fragL4_dif =  SSC.add_diffusion(fragL4, T, B, G, M)
    fragL44_dif =  SSC.add_diffusion(fragL44, T, B, G, M)
    fragL100_dif =  SSC.add_diffusion(fragL100, T, B, G, M)

    print "--BD standard deviation + diffusion error--"
    print msg.format(1000, np.std(fragL1_dif))
    print msg.format(4000, np.std(fragL4_dif))
    print msg.format(44500, np.std(fragL44_dif))
    print msg.format(100000, np.std(fragL100_dif))

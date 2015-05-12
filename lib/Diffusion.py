#!/usr/bin/env python

# import
## batteries
import sys,os
import cPickle as pickle 
from functools import partial
## 3rd party
import parmap
import scipy.stats as stats
import numpy as np
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from FragGC import Frag_multiKDE
import SIPSimCython as SSC
import Utils


# functions
def make_kde(taxon_name, kde, n, bw_method, **kwargs):
    if kde is None:
        return (taxon_name, None)

    BD_wDiff = SSC.add_diffusion(kde.resample(size=n),
                                 **kwargs)

    kdeBD = stats.gaussian_kde(BD_wDiff,
                               bw_method=bw_method)

    return (taxon_name, kdeBD)
    

def main(args):    
    """
    load kde object
    foreach kde (parallel):
      foreach iter:
        draw from kde and calculate diffusion
        append to array
      calculate BD from kde
      create kde from BD array
      return kde object
    """
    try:
        args['--bw'] =float(args['--bw'])
    except TypeError:
        pass 

    kde2d = Utils.load_kde(args['<fragment_kde>'])

    pfunc = partial(make_kde, 
                    n = int(args['-n']),
                    T = float(args['-T']),
                    B = float(args['-B']),
                    G = float(args['-G']),
                    bw_method=args['--bw'])
    
    KDE_BD = parmap.starmap(pfunc, kde2d.items(),
                            processes = int(args['--np']),
                            chunksize = int(args['--cs']),
                            parallel = not args['--debug'])
    

    pickle.dump({taxon:KDE for taxon,KDE in KDE_BD}, 
                sys.stdout)
    


        
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

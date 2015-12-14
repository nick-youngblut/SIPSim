#!/usr/bin/env python

# import
## batteries
import sys,os
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

import SIPSimCython as SSC
import Utils

# functions
def BD2DBL_index(BD, step=0.1):
    """Based on provided gradient params, make a dict that
    relates DNA fragment GC to the GC span equilent for the DBL.
    For instance, a GC of 50 (%) may make a DBL spanning the gradient
    location equivalent of 30-70% GC, and thus the 50% GC fragments
    in the DBL may end up anywhere in 30-70% GC.

    Returns
    -------
    dict : (BD : DBL)
           DBL = uniform distribution function that returns new BD value
    """
#    GCs = np.arange(0,100,step)  # GC range
#    BDs = GC2BD(GCs)
    BDs = np.arange(BD_min, BD_max, BD_step)
    angle_tube_pos = BD2angleTubePos(BDs)
    angle_tube_vol = BD2angledTubeVol(BDs)
    vert_tube_pos = angledTubeVol2vertTubeHeight(angle_tube_vols)    
    VTP2BD_func = vertTubePos_fit(vert_tube_pos, BDs)
    DBL = get_DBL(angle_tube_pos, BDs, VTP2BD_func)
    #DBL = DBL_BD2GC(DBL)
  

def BD2angleTubePos(BD):
    """Convert BD to angled tube position
    """
    pass

def BD2angledTubeVol(BD):
    """Convert BD to volume of gradient where the BD is >= to the provided BD.
    """
    pass

def angledTubeVol2vertTubeHeight(v):
    """Convert angled tube volume (see BD2angledTubeVol) to height in the
    vertical tube.
    """
    pass

#def make_BD_pos_index(BDs, angle_tube_pos, vert_tube_pos):
#    """Make index for each BD (BD : (angle|vert) : [min_pos/max_pos | pos])
#    """
#    pass


def vertTubePos_BD_fit(vert_tube_pos, BDs):
    """Making a function: vert_tube_pos ~ BD.
    Using curve fit to define function.
    This function will be used to convert angle_tube_pos to 
    vert_tube_BD
    """
    pass

def get_DBL(angle_tube_pos, BDs, VTP2BD_func):
    """Converting angle_tube_pos (min/max) to vertical tube BD
    min/max (span of DBL), converting span to uniform distribution
    function (min/max of DBL), then making index: (BD : uniform_dist_func)
    """
    pass


def KDE_with_DBL(BD_KDE, DBL_index, n, frac):
    """Sample <frac> from BD_KDE and convert values to DBL BD values
    and sample 1 - <frac> from BD_KDE; combine all BD values and make a KDE
    
    Returns
    -------
    KDE
    """
    # input unpacking a type checking                          
    try:
        taxon_name,kde = x   
    except ValueError:
        msg = '"x" must be (taxon_name, kde)'
        raise ValueError, msg   
    try:
        bw_method = float(bw_method)                                            |
    except (TypeError, ValueError) as e:
        pass                                                                    |
    try:
        bw_method = bw_method.lower()                                           |
    except AttributeError:
        pass


    # status
    msg = 'Processing: {}\n'
    sys.stderr.write(msg.format(taxon_name))                                    |
      
    # if no KDE
    if kde is None: 
        return (taxon_name, None)
        
    # sampling fraction for DBL

    # sampling 1-DBL-fraction for non-DBL
    
    # making new KDE from samplings
    #kdeBD = stats.gaussian_kde(BD_vals, bw_method=bw_method)


def main(args):    
    """Main function for adding G+C value error due to DBL
    Parameters:
    args : dict
        See ``DBL`` subcommand
    """
    # making BD2DBL index
    

    # loading fragment KDEs of each genome
    kde2d = Utils.load_kde(args['<fragment_kde>'])

    # for each genome KDE, making new KDE with DBL 'smearing'
    ## for '-n' selecting fraction to be in DBL, other is in standard
    ### for DBL fraction, sampling from KDE, get new BD values
    ### for non-DBL fraction, sampling from KDE, append to DBL BD values
    ## make a new KDE of BD values

#    pfunc = partial(KDE_with_DBL, 
#                    n = int(args['-n']),
#                    T = float(args['-T']),
#                    B = float(args['-B']),
#                    G = float(args['-G']),
#                    bw_method=args['--bw'])
    
    # difussion calc in parallel
#    pool = ProcessingPool(nodes=int(args['--np']))
#    if args['--debug']:
#        KDE_BD = map(pfunc, kde2d.items())
#    else:
#        KDE_BD = pool.map(pfunc, kde2d.items())

    # pickling output
#    dill.dump({taxon:KDE for taxon,KDE in KDE_BD}, sys.stdout)
    


        
if __name__ == '__main__':
    pass

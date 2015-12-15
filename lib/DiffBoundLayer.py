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


# functions for trig with angle in degrees
cos_d = lambda d : np.cos(np.deg2rad(d))
sin_d = lambda d : np.sin(np.deg2rad(d))
asin_d = lambda d : np.arcsin(np.deg2rad(d))
acos_d = lambda d : np.arccos(np.deg2rad(d))


def calc_tube_angle(r_min, r_max, tube_height):
    """Convert min/max distance from axis of rotation
    (+ tube height) to angle of tube (units = degrees).
    
    Parameters
    ----------
    r_min : float
            Min distance to axis of rotation (cm)
    r_max : float
            Max distance to axis of rotation (cm)    
    tube_height : float
            Height of cfg tube (cm)
    """
    x = r_max - r_min
    hyp = tube_height
    rotor_angle = np.rad2deg(np.arcsin(x / hyp))
    return rotor_angle


def calc_isoconc_point(r_min, r_max):
    """Convert min/max distance from axis of rotation
    to the isoconcentration point of the gradient
    (units = cm).
    
    Parameters
    ----------
    r_min : float
            min distance to axis of rotation (cm)
    r_max : float
            max distance to axis of rotation (cm)    
    """
    I = np.sqrt((r_min**2 + r_min * r_max + r_max**2)/3)
    if not isinstance(I, float):
        msg = 'isoconcentration point calc error: {} is not a float.'
        raise TypeError, msg.format(I)
    return I


def BD2distFromAxis(BD, D, B, w, I):
    """Convert buoyant density (BD) to distance from axis of 
    rotation (units = cm).
    
    Parameters
    ----------
    BD : float
        Buoyant density of particle
    D : float
        Average density of the gradient
    B : float
        Beta coefficient (Beta^o) 
    w : float
        Angular velocity (omega^2)
    I : float
        Isoconcentration point (cm)
    """
    X = np.sqrt(((BD-D)*2*B/w) + I**2)
    if not isinstance(X, float):
        msg = 'distFromAxis calc error: {} is not a float.'
        raise TypeError, msg.format(X)
    return X
    

def axisDist2angledTubePos(x, tube_radius, r_max, A):
    """Convert distance from axis of rotation to angled tube position
    
    Parameters
    ----------
    x : float
        Distance from axis of rotation for a particle (cm)
    tube_radius : float
        Radius of cfg tube (cm)
    r_max : float
        Max cfg tube distance from axis of rotation (cm)
    A : angle of tube to axis of rotation (degrees)
    """
    if(x >= r_max - (tube_radius * cos_d(A)) - tube_radius):
        # band in rounded bottom of cfg tube
        d = x - (r_max - tube_radius)
        a = A - asin_d(d / tube_radius)
        LowH = tube_radius - tube_radius * cos_d(a)
        #print LowH
    else:
        # band in cylinder of cfg tube
        d = r_max - (tube_radius * cos_d(A)) - tube_radius - x
        h_c = d/sin_d(A)
        LowH = tube_radius + h_c
        # print LowH

    if(x > r_max - (tube_radius - tube_radius * cos_d(A))):
        # Equation for finding the upper band
        d = x - (r_max - tube_radius)
        a = A - (180 - asin_d(d/tube_radius))
        HighH = tube_radius - tube_radius * cos_d(a)
        #print HighH
    else:
        # This band will be in the cylinder part
        d = r_max - (tube_radius - tube_radius * cos_d(A)) - x
        h_c = d/sin_d(A)
        HighH = tube_radius + h_c
        #print(HighH)     
    return(LowH, HighH)


def axisDist2angledTubeVol(axisDist):
    """Convert distance from axis of rotation to volume of gradient
    where the BD is >= to the provided BD.
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
    """Making a continuous function: vert_tube_pos ~ BD.
    Using curve fit to define function.
    This function will be used to convert angle_tube_pos to 
    vert_tube_BD
    """

    # scipy curve_fit
    ## http://www2.mpia-hd.mpg.de/~robitaille/PY4SCI_SS_2014/_static/15.%20Fitting%20models%20to%20data.html
    ## lm_func = lambda x,a,b : a*x + b
    ## pm2_func = lambda x,a,b,c : a*x**2 + b*x + c
    ## pm3_func = lambda x,a,b,c,d : a*x**3 + b*x**2 + c*x + d
    # np.polyfit?

    pass


def get_DBL(angle_tube_pos, BDs, VTP2BD_func):
    """Converting angle_tube_pos (min/max) to vertical tube BD
    min/max (span of DBL), converting span to uniform distribution
    function (min/max of DBL), then making index: (BD : uniform_dist_func)
    """
    # np.random.uniform
    pass


def BD2DBL_index(r_min, r_max, D, B, w,
                 tube_diam, tube_height,
                 BD_min=1.67, BD_max=1.78, BD_step=0.001):
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
    # vectorizing functions
    BD2distFromAxisV = np.vectorize(BD2distFromAxis)
    axisDist2angledTubePosV = np.vectorize(axisDist2angledTubePos)
    axisDist2angledTubeVolV = np.vectorize(axisDist2angledTubeVol)

    # calculate gradient/cfg_tube params
    tube_radius = tube_diam / 2
    I = calc_isoconc_point(r_min, r_max)    
    A = calc_tube_angle(r_min, r_max, tube_height)

    # BD value range
    BDs = np.arange(BD_min, BD_max, BD_step)

    # convert to distance from axis of rotation
    axisDists = BD2distFromAxisV(BDs, D, B, w, I)
        
    # angle tube position/volume info
    angle_tube_pos = axisDist2angledTubePosV(axisDists, 
                                             tube_radius, 
                                             r_min, A)
    #-- debug
    return None
    #print angle_tube_pos
    
    angle_tube_vol = axisDist2angledTubeVolV(axisDists)
    
    # determine the continuous function: BD_vertTube ~ vert_tube_height
    vert_tube_pos = angledTubeVol2vertTubeHeight(angle_tube_vols)    
    VTP2BD_func = vertTubePos_fit(vert_tube_pos, BDs)

    # convert angle tube position info to BD range of DBL
    DBL_index = get_DBL(angle_tube_pos, BDs, VTP2BD_func)
    return DBL_index



def KDE_with_DBL(BD_kde, DBL_index, n, frac_abs, bw_method):
    """Sample <frac> from BD_KDE and convert values to DBL BD values
    and sample 1 - <frac> from BD_KDE; combine all BD values and make a KDE
    
    Returns
    -------
    KDE of BD values
    """
    # input unpacking a type checking                          
    try:
        taxon_name,kde = BD_kde
    except ValueError:
        msg = '"BD_kde" must be (taxon_name, kde)'
        raise ValueError, msg   
    try:
        bw_method = float(bw_method)
    except (TypeError, ValueError) as e:
        pass
    try:
        bw_method = bw_method.lower()
    except AttributeError:
        pass
    msg = '"frac_abs" must be a fraction (0-1)'
    assert 0 <= frac_abs <= 1, msg
    msg = 'DBL_index must be a vectorized function'
    assert isinstance(DBL_index, np.lib.function_base.vectorize), msg
    

    # status
    msg = 'Processing: {}\n'
    sys.stderr.write(msg.format(taxon_name)) 
      
    # if no KDE
    if kde is None:
        return (taxon_name, None)

    # sampling fraction for DBL
    n_DBL = int(n * frac_abs)
    roundV = np.vectorize(round)
    DBL_frags = roundV(kde.resample(size=n_DBL), 3)    
    DBL = DBL_index(DBL_frags)    

    # sampling (1 - DBL_fraction) for non-DBL
    n_nonDBL = int(n - n_DBL)
    x = kde.resample(size=n_nonDBL)
    print DBL.shape
    print x.shape
    #exit()
    BD_new = np.concatenate((DBL, kde.resample(size=n_nonDBL)), axis=1)

    print BD_new; exit()
 
    # making new KDE from samplings
#    kdeBD = stats.gaussian_kde(np.concatenate(
#                                               kde.resample(size=n_nonDBL))),
#                               bw_method=bw_method)
#    print kdeBD; exit()


def main(args):    
    """Main function for adding G+C value error due to DBL
    Parameters:
    args : dict
        See ``DBL`` subcommand
    """

    # making BD2DBL index
    DBL_index = BD2DBL_index(r_min = float(args['--r_min']),
                             r_max = float(args['--r_max']),
                             D = float(args['-D']),
                             B = float(args['-B']),
                             w = float(args['-w']),
                             tube_diam = float(args['--tube_diam']),
                             tube_height = float(args['--tube_height']),
                             BD_min = float(args['--BD_min']),
                             BD_max = float(args['--BD_max']))

    #--debug--#
    ## pseudo DBL_index
    BDs = np.arange(1.64, 1.78, 0.001)
    DBL_index = {}
    for x, y, z in zip(BDs, BDs, BDs):
        x = round(x, 3)
        DBL_index[x] = partial(np.random.uniform,
                               low = y,
                               high = z + 0.001,
                               size=1)
        
    #for x in BDs:
    #    x = round(x, 3)
    #    print DBL_index[x](**{'size' : 2})


    # loading fragment KDEs of each genome
    kde2d = Utils.load_kde(args['<fragment_kde>'])

    # for each genome KDE, making new KDE with DBL 'smearing'
    ## for '-n' selecting fraction to be in DBL, other is in standard
    ### for DBL fraction, sampling from KDE, get new BD values
    ### for non-DBL fraction, sampling from KDE, append to DBL BD values
    ## make a new KDE of BD values
    DBL_indexV = np.vectorize(lambda x : DBL_index[x]())
    pfunc = partial(KDE_with_DBL, 
                    DBL_index = DBL_indexV,
                    n = int(args['-n']),
                    frac_abs = float(args['--frac_abs']),
                    bw_method = args['--bw'])
    
    # difussion calc in parallel
    pool = ProcessingPool(nodes=int(args['--np']))
    if args['--debug']:
        KDE_BD = map(pfunc, kde2d.items())
    else:
        KDE_BD = pool.map(pfunc, kde2d.items())

    # pickling output
#    dill.dump({taxon:KDE for taxon,KDE in KDE_BD}, sys.stdout)
    


        
if __name__ == '__main__':
    pass

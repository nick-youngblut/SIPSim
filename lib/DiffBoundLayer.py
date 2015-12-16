#!/usr/bin/env python

# import
## batteries
import sys,os
from functools import partial
## 3rd party
import scipy.stats as stats
from scipy.integrate import quad
from scipy.optimize import brentq
import numpy as np
import dill
from pathos.multiprocessing import ProcessingPool
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
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

    Returns
    -------
    float : distance from axis (cm)
            Note: 
    """
    X2 = ((BD - D) * 2*B/w) + I**2
    if X2 < 0:
        return np.nan
    X = np.sqrt(X2)
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

    Returns
    -------
    tuple of lowest & highest positions in tube vertical height that the 
    band extends to. Note: nan returned if x = nan
    """
    if np.isnan(x):
        return (x, x)

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


def _SphVol(t, r, p2, R12):
    # helper function for axisDist2angledTubeVol
    v1 = t*((2*r)-t)/2
    v2 = 2*np.pi*((p2-t)/R12)
    v3 = np.sin(2*np.pi*((p2-t)/R12))
    return v1 * (v2 - v3)
    
def _CylWedVol(t, r, b, h):
    # helper function for axisDist2angledTubeVol
    return 2*(h*(t-r+b)/ b) * np.sqrt(r**2-t**2)

def axisDist2angledTubeVol(x, r, D, A):
    """Convert distance from axis of rotation to volume of gradient
    where the BD is >= to the provided BD.
    
    Parameters
    ----------
    x : float
        distance from axis of rotation (cm)    
    r : float
        cfg tube radius (cm)
    D : float
        max distance from axis of rotation (cm)
    A : float
        cdf tube angle in rotor (degrees)
    
    Returns
    -------
    volume (ml) occupied by gradient heavier or as heavy as at that point.
    Note: nan returned if x = nan
    """
    # return nan if nan provided
    if np.isnan(x):
        return x

    a = np.deg2rad(A)
    p1 = r-(r*np.cos(a))
    p2 = r+(r*np.cos(a))
    R12 = p2-p1
    d = D-x
    D1 = D-p1
    D2 = D-p2        
    msg = ' '.join(['WARNING: Cannot compute angled tube volume!', 
                    '{} is beyond end of tube.\n'])

    if x < D2:
        h1 = (D2-x)/np.sin(a)
        h2 = (D1-x)/np.sin(a)
        volume1 = (2/3.0)*np.pi*r**3
        volume2 = (0.5)*np.pi*r**2*(h1+h2)
        volume = volume1+volume2
    elif D1 >= x >= D2:
        volume1 = (1/3.0)*np.pi*p1**2*(3*r-p1)
        volume2 = quad(_SphVol, p1, d, args=(r, p2, R12))
        b = (d-p1)/np.cos(a)
        h = b/np.tan(a)
        volume3 = quad(_CylWedVol, r-b, r, args=(r, b, h))
        volume = volume1+volume2[0]+volume3[0]
    elif D >= x > D1:
        volume = (1/3.0)*np.pi*d**2*(3*r-d)
    elif x > D:
        sys.stderr.write(msg.format(x))
        volume = np.nan
        #raise ValueError, msg.format(x)
    else:
        sys.stderr.write(msg.format(x))
        volume = np.nan
        #raise ValueError, msg.format(x)
        
    return volume


def _cylVol2height(v, r): 
    # v = volume (ml)
    # r = tube radius (cm)
    h = v / (np.pi * r**2)
    return h

def _sphereCapVol2height(v, r):
    # v = volume (ml)
    # r = tube radius (cm)
    # height = h**3 - 3*r*h**2 + (3v / pi) = 0
    f = lambda x : x**3 - 3*r*x**2 + 3*v/np.pi
    try:
        root = brentq(f, 0, r*2, maxiter=1000)
    except ValueError:
        msg = 'WARNING: not roots for volume {}\n'
        sys.stderr.write(msg.format(v))
        root = np.nan
    return(root)


def tubeVol2vertTubeHeight(v, r):
    """Convert angled tube volume (see axisDist2angledTubeVol) to height in the
    vertical tube.
    
    Parameters
    ----------
    v : float
        Volume (ml)
    r : float
        Tube radius (cm)
    """
    sphere_cap_vol = (4/3 * np.pi * r**3)/2
    
    if v <= sphere_cap_vol:
        # height does not extend to cylinder
        h = _sphereCapVol2height(v, r)
    else:
        # height = sphere_cap + cylinder
        sphere_cap_height = _sphereCapVol2height(sphere_cap_vol, r)
        h =  sphere_cap_height + _cylVol2height(v - sphere_cap_vol, r)

    return(h)


def vertTubePos_BD_fit(BDs, vert_tube_pos, deg=3):
    """Making a continuous function: BD ~ vert_tube_pos.
    This function will be used to convert angle_tube_pos to 
    the BD in the vertical tube gradient. Using polynomial curve fit to 
    define function.
    
    Parameters
    ----------
    BDs : 1d array
         Buoyant density values
    vert_tube_pos : 1d array
         Vertical tube height corresponding to BD at same index
    deg : int
         Degree of the fitting polynomial
    
    Returns
    -------
    numpy.poly1d instance
    """
    deg = int(deg)
    # must have equal length arrays
    msg = 'Number of vertical tube positions != number of BD values'
    assert len(vert_tube_pos) == len(BDs), msg
    
    # trimming NA values
    vert_tube_pos = vert_tube_pos[~np.isnan(vert_tube_pos)]
    BDs = BDs[~np.isnan(vert_tube_pos)]

    # fitting variables (polynomial fitting)
    fit = np.polyfit(vert_tube_pos, BDs, deg=deg)
    return  np.poly1d(fit)


def make_DBL_index(angle_tube_pos, BDs, VTP2BD_func):
    """Converting angle_tube_pos (min/max) to vertical tube BD
    min/max (span of DBL), converting span to uniform distribution
    function (min/max of DBL), then making index: (BD : uniform_dist_func)
    
    Parameters
    ----------
    angle_tube_pos : tuple of 1d arrays
         ([array_of_lowest_angled_tube_pos],[array_of_highest_angled_tube_pos]).
    BDs : 1d array
         Buoyant density values corresponding to angle_tube_pos at each index.
    VTP2BD_func : callable
         Function that accepts a 1d array (angle_tube_pos) and returns BD values
         for same tube position in vertical tube gradient.
    """
    # must have equal length arrays
    msg = 'Number of angle tube positions != number of BD values'
    assert len(angle_tube_pos) == 2, 'Formatting error: angled tube positions'
    assert len(angle_tube_pos[0]) == len(BDs), msg
    assert len(angle_tube_pos[1]) == len(BDs), msg

    # converting angle_tube_pos values to BD values
    low_pos_BD = VTP2BD_func(angle_tube_pos[0])
    high_pos_BD = VTP2BD_func(angle_tube_pos[1])
    
    # making a dict of BD : np.random.uniform(low_pos_BD, high_pos_BD)
    DBL_index = {}
    for i in xrange(len(low_pos_BD)):
        BD = round(BDs[i], 3)
        if np.isnan(low_pos_BD[i]) or np.isnan(high_pos_BD[i]):
            pass
        else:
            DBL_index[BD] = partial(np.random.uniform,
                                    low = low_pos_BD[i],
                                    high = high_pos_BD[i])
    return DBL_index


def BD2DBL_index(r_min, r_max, D, B, w,
                 tube_diam, tube_height,
                 BD_min=1.67, BD_max=1.78, BD_step=0.001):
    """Based on provided gradient params, make a dict that
    relates DNA fragment GC to the GC span equilent for the DBL.
    For instance, a GC of 50 (%) may make a DBL spanning the gradient
    location equivalent of 30-70% GC, and thus the 50% GC fragments
    in the DBL may end up anywhere in 30-70% GC. 

    Parameters
    ----------

    Returns
    -------
    dict : (BD : DBL)
           DBL = uniform distribution function that returns new BD value
    """
    # vectorizing functions
    BD2distFromAxisV = np.vectorize(BD2distFromAxis)
    axisDist2angledTubePosV = np.vectorize(axisDist2angledTubePos)
    axisDist2angledTubeVolV = np.vectorize(axisDist2angledTubeVol)
    tubeVol2vertTubeHeightV = np.vectorize(tubeVol2vertTubeHeight)

    # calculate gradient/cfg_tube params
    tube_radius = tube_diam / 2
    I = calc_isoconc_point(r_min, r_max)    
    A = calc_tube_angle(r_min, r_max, tube_height)

    # BD value range
    BDs = np.arange(BD_min, BD_max, BD_step)

    # convert to distance from axis of rotation
    axisDists = BD2distFromAxisV(BDs, D, B, w, I)
        
    # angle tube position/volume info
    angle_tube_pos = axisDist2angledTubePosV(axisDists, tube_radius, r_max, A)
    angle_tube_vol = axisDist2angledTubeVolV(axisDists, tube_radius, r_max, A)
    
    # determine the continuous function: BD_vertTube ~ vert_tube_height
    vert_tube_pos = tubeVol2vertTubeHeightV(angle_tube_vol, tube_radius)    
    VTP2BD_func = vertTubePos_BD_fit(BDs, vert_tube_pos)

    # convert angle tube position info to BD range of DBL
    DBL_index = make_DBL_index(angle_tube_pos, BDs, VTP2BD_func)
    return DBL_index


def _fake_DBL_index(BD_min, BD_max, BD_step=0.001):
    # pseudo DBL_index
    BDs = np.arange(BD_min, BD_max, BD_step)
    DBL_index = {}
    for x, y, z in zip(BDs, BDs, BDs):
        x = round(x, 3)
        DBL_index[x] = partial(np.random.uniform,
                               low = y,
                               high = z + 0.001,
                               size=1)
    return DBL_index


def KDE_with_DBL(BD_kde, DBL_index, n, frac_abs, bw_method):
    """Sample <frac> from BD_KDE and convert values to DBL BD values
    and sample 1 - <frac> from BD_KDE; combine all BD values and make a KDE
    
    Parameters
    ----------
    BD_kde : list
         [taxon_name, scipy_gausiuan_kde_object]
    DBL_index : dict
         {BD : distribution_function_that_returns_a_value}
    n : int
         Total number of subsampling from KDE + DBL
    frac_abs : float 
         Fraction of DNA fragments in DBL.  range = (0-1)
    bw_method : see stats.gaussian_kde

    Returns
    -------
    stats.gaussian_kde object describing distribution of fragment BD values 
    (DBL 'smearing' included)
    """
    n = int(n)
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
    DBL_frags = np.round(kde.resample(size=n_DBL), 3)    
    DBL_frags = DBL_index(DBL_frags)    

    # sampling (1 - DBL_fraction) for non-DBL
    n_nonDBL = int(n - n_DBL)
 
    # making new KDE from samplings
    kdeBD = stats.gaussian_kde(np.concatenate((DBL_frags, 
                                               kde.resample(size=n_nonDBL)), 
                                              axis=1)[0,],
                               bw_method=bw_method)
    return (taxon_name, kdeBD)


def main(args):    
    """Main function for adding DBL 'smearning' to fragment distribution in
    the CsCl gradient.

    Parameters
    ----------
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
    #DBL_index = _fake_DBL_index(BD_min = float(args['--BD_min']),
    #                            BD_max = float(args['--BD_max']))
        
    # loading fragment KDEs of each genome
    kde2d = Utils.load_kde(args['<fragment_kde>'])

    # for each genome KDE, making new KDE with DBL 'smearing'
    def DBL_indexF(x):
        msg = 'WARNING: BD value "{}" not found in DBL\n'
        try:
            y = DBL_index[x]()
        except KeyError:
            sys.stderr.write(msg.format(x))
            x = np.random.choice(DBL_index.keys())
            y = DBL_index[x]()
        return y
    
    DBL_indexFV = np.vectorize(DBL_indexF)                              
    pfunc = partial(KDE_with_DBL, 
                    DBL_index = DBL_indexFV,
                    n = int(args['-n']),
                    frac_abs = float(args['--frac_abs']),
                    bw_method = args['--bw'])
    
    # DBL KDE calc in parallel (per taxon)
    pool = ProcessingPool(nodes=int(args['--np']))
    if args['--debug']:
        KDE_BD = map(pfunc, kde2d.items())
    else:
        KDE_BD = pool.map(pfunc, kde2d.items())

    # pickling output
    dill.dump({taxon:KDE for taxon,KDE in KDE_BD}, sys.stdout)    

        
if __name__ == '__main__':
    pass

# import
## batteries
import os, sys
import re
from functools import partial
import types
import logging
from collections import defaultdict
from random import shuffle
from StringIO import StringIO
## 3rd party
import numpy as np
import pandas as pd
import mixture
import scipy.stats as stats
import dill
from pathos.multiprocessing import ProcessingPool
## application
import Utils
import SIPSimCython
from CommTable import CommTable
from SIPSimCython import add_incorp
from Config import Config

# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def main(args):
    """Main function for performing isotope incorporation simulation.
    """
    # loading input
    ## kde_bd
    KDE_BD = Utils.load_kde(args['<BD_KDE>'])
    ## comm (optional)
    if args['--comm'] is not None:
        comm = CommTable.from_csv(args['--comm'], sep='\t')
        # combining kde and comm
        _add_comm_to_kde(KDE_BD, comm)
    else:
        comm = None

    
    # loading the config file
    config = Config.load_config(args['<config_file>'],
                                phylo=args['--phylo'])
     

    # creating kde of BD distributions with BD shift from isotope
    KDE_BD_iso = dict()
    for libID in config.keys():        
        # making a list of taxa that can incorporate 
        # TODO: abundance cutoff: taxa must have abundance > threshold to incorporate
        taxa_incorp_list = _taxon_incorp_list(libID, config, KDE_BD)

        # TODO: abundance weighting with less incorp for less taxa

        # setting params for parallelized function
        pfunc = partial(_make_kde, 
                        config = config,
                        libID = libID,
                        n = args['-n'],
                        taxa_incorp_list = taxa_incorp_list,
                        isotope = args['--isotope'],
                        bw_method = args['--bw'])                    

        # parallel by taxon
        pool = ProcessingPool(nodes=int(args['--np']))
        if args['--debug']:
            tmp = map(pfunc, KDE_BD.items())
        else:
            tmp = pool.map(pfunc, KDE_BD.items())

        KDE_BD_iso[libID] = {taxon:KDE for taxon,KDE in tmp}
    

    # writing pickled BD-KDE with BD shift from isotope incorp
    dill.dump(KDE_BD_iso, sys.stdout)
        

def _make_kde(x, libID, config, taxa_incorp_list, 
             isotope='13C', n=10000, bw_method=None): 
    """Making new KDE of BD value distribution which includes
    BD shift due to isotope incorporation. 
    Args:
    x -- [taxon_name, dict -- {kde:abundance}]
    libID -- str; library ID
    config -- config object
    taxa_incorp_list -- iterable; taxa that can incorporate isotope
    isotope -- str; isotope that is incorporated
    n -- number of Monte Carlo samples to use for estimating
         BD+isotope_BD distribution
    bw_method -- bandwidth scalar or function passed to
                 scipy.stats.gaussian_kde().
    Return:
    (taxon_name, KDE*)
       * Note: KDE object may be None    
    """
    taxon_name,x = x

    # status
    sys.stderr.write('Processing: {}\n'.format(taxon_name))

    # unpack
    n = int(n)
    kde = x['kde']
    if kde is None:
        return (taxon_name, None)

    # can taxon incorporate any isotope?
    ## if not, return 'raw' BD KDE 
    if taxon_name not in taxa_incorp_list:
        return (taxon_name, kde)

    # taxon abundance for library
    try:
        taxon_abund = x['abundances'][libID]
    except KeyError:
        tmp = re.sub('[Ll]ib(rary)* *','', libID)
        try:
            taxon_abund = x['abundances'][tmp]
        except KeyError:
            taxon_abund = None

    # max incorporation for isotope
    maxIsotopeBD = isotopeMaxBD(isotope)

    # making a mixture model object lib:taxon
    ## mix_model => distribution of % incorporation for taxon population
    mix_model = make_incorp_model(taxon_name, libID, config)

    # making KDE of BD + BD_isotope_incorp
    kdeBD = stats.gaussian_kde( 
        add_incorp(kde.resample(n)[0], 
                   mix_model, 
                   maxIsotopeBD),
        bw_method=bw_method)

    return (taxon_name, kdeBD)


def _add_comm_to_kde(KDE_BD, comm):
    """Adding comm data for each taxon to each KDE.
    Args:
    Return:
    KDE_BD -- {taxon_name:{kde|abundances}}
    """
    libIDs = comm.get_unique_libIDs()
    for taxon_name,kde in KDE_BD.items():
        d = {'kde':kde, 'abundances':{}}
        for libID in libIDs:
            abund = comm.get_taxonAbund(taxon_name, libID=libID)
            d['abundances'][libID] = abund[0]
        KDE_BD[taxon_name] = d    


def _taxon_incorp_list(libID, config, KDE_BD):
    # perc incorp from config
    try:
        max_perc_taxa_incorp = config.get_max_perc_taxa_incorp(libID) 
    except KeyError:
        max_perc_taxa_incorp = 100.0
    max_perc_taxa_incorp /= 100.0

    # randomized list of taxa
    taxon_names = KDE_BD.keys()
    shuffle(taxon_names)

    # set of taxa that incorporate any isotope
    n_incorp = int(round(len(taxon_names) * max_perc_taxa_incorp, 0))
    return taxon_names[:n_incorp]

    

def make_incorp_model(taxon_name, libID, config):
    """Setting isotope incorporation based on the interPopDist
    function for each intraPop parameter.
    Args:
    taxon_name -- taxon name string
    libID -- library ID string
    config -- config object
    Returns:
    mixture model object (mixture class)
    """
    psblDists = {'normal' : mixture.NormalDistribution,
                 'uniform' : mixture.UniformDistribution}

    # creating individual distribution functions
    libSect = config.get_libSection(libID)
    intraPopDist_IDs = []
    intraPopDist_funcs = []
    weights = []
    for (intraPopDistID,intraPopDist) in config.iter_configSections(libSect):    

        # name of standard distribution
        try:
            distID = intraPopDist['distribution']
        except KeyError:
            msg = 'Cannot find "distribution" key for "{}"'
            raise KeyError, msg.format(intraPopDistID)
        intraPopDist_IDs.append(distID)
        # getting mixture model weight for this distribution
        weight = float(intraPopDist.get('weight', 0))
        weights.append(weight)
        # selecting intra-pop param values from inter-pop dists
        params = dict()        
        for (paramID,param) in config.iter_configSections(intraPopDist):
            params[paramID] = _select_intrapop_param_value(param, 
                                                           taxon_name)

        # checking start-end parameters (if present)
        _start_lt_end(params)

        # making intra-pop dist function (a standard distribution)
        try:
            dist_func = psblDists[distID](**params)
        except KeyError:
            msg = 'Distribution "{}" not supported'
            raise KeyError, msg.format(distID)
        intraPopDist_funcs.append(dist_func)
    
    
    # making sure weights add up to 1
    weights = Config._fill_in_weights(weights)
    assert len(weights) == len(intraPopDist_IDs), \
        'number_of_distributions != number_of_weights'
    assert len(intraPopDist_IDs) == len(intraPopDist_funcs), \
        'number_of_distributions != number_of_dist_functions'

    # making incorporation mixture model
    return mixture.MixtureModel(len(intraPopDist_IDs),
                                weights,
                                intraPopDist_funcs)
                            


def _select_intrapop_param_value(interPopDist, taxon_name, maxtries=1000):
    """Selecting the intra-population parameter value
    based on the inter-population distribution function.
    Values are % isotope incorporation, so acceptable
    range is 0-100 (will try 'maxtries' times to select value in range).
    Args:
    interPopDist -- {'interPopDist':{'function':interPopdist_function}}
    taxon_name -- name of taxon
    maxtries -- number of tries to get a parameter values >0
    Returns:
    float -- intra-pop param value
    """
    # getting inter-pop dist function
    try:
        interPopDist_func = interPopDist['interPopDist']['function']
    except KeyError:
        raise KeyError, 'Cannot find inter-pop dist function'


    # sampling from function to get parameter for intra-pop distribution
    tries = 0
    while True:
        tries += 1
        try:
            # if Brownian Motion evolution
            paramVal = interPopDist_func.sample(taxon_name)
        except TypeError:
            paramVal = interPopDist_func.sample()
            try:
                paramVal = paramVal[0]
            except TypeError:
                pass
        # values must be >= 0 and <= 100 to end loop
        if paramVal >= 0 and paramVal <= 100:
            break
        # exceeding maxtries?
        if tries >= maxtries:
            err = 'Taxon: {}'.format(taxon_name)
            msg = 'Exceeded maxtries to get parameter in range: 0-100'  
            sys.exit(': '.join([err, msg]))

    return paramVal


    
def _start_lt_end(params):
    """Check that 'start' param is < 'end' param
    if both params are found in the provided dict.
    Args:
    params -- {param_ID:param_value}
    Returns:
    in-place edit of params
    """
    if ('start' in params) & ('end' in params):
        try:
            startVal = float(params['start'])
            endVal = float(params['end'])
        except TypeError:
            return None            
            
        if startVal > endVal:
            params['start'] = endVal
            params['end'] = startVal
        elif startVal == endVal:
            # start-end cannot ==
            if endVal >= 100:
                params['start'] = startVal - 1e-10
            else:
                params['end'] = endVal + 1e-10


def isotopeMaxBD(isotope):
    """Setting the theoretical max BD shift of an isotope (if 100% incorporation).
    Args:
    isotope -- str; name of isotope
    Return:
    float -- max BD value
    """
    psblIsotopes = {'13C' : 0.036,
                    '15N' : 0.016}
    try:
        return psblIsotopes[isotope.upper()]
    except KeyError:
        raise KeyError('Isotope "{}" not supported.'.format(isotope))

        

    

# import
## batteries
import os, sys
import re
from functools import partial
import types
import logging
import tempfile
import shutil
import glob
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
from CommTable import CommTable
from IsoIncorpCython import add_incorp
from Config import Config

# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def main(args):
    """Main function for performing isotope incorporation simulation.
    
    Parameters
    ----------
    args : dict
        See ``isotope_incorp`` subcommand.
    """
    # loading input
    ## config file
    config = Config.load_config(args['<config_file>'],
                                phylo=args['--phylo'])

    ## comm (optional)
    if args['--comm'] is not None:
        comm = CommTable.from_csv(args['--comm'], sep='\t')
    else:
        comm = None

    ## taxa list (optional
    if args['--taxa'] is not None:
        taxa_incorp_list = _load_taxa_incorp_list(args['--taxa'], config)
    else:
        taxa_incorp_list = {k:[] for k in config.keys()}
            
    ## kde_bd
    sys.stderr.write('Loading KDE object...\n')
    KDEs = Utils.load_kde(args['<BD_KDE>'])

    # which depth of KDE object
    KDE_type = Utils.KDE_type(KDEs)

    # adding comm info to KDEs 
    if KDE_type < 4:
        _add_comm_to_kde(KDEs, comm)
        
    # temporary directory for BD shift stats
    tmpdirpath = tempfile.mkdtemp()    

    # creating kde of BD distributions with BD shift from isotope
    KDEs_iso = dict()
    for libID in config.keys():        
        sys.stderr.write('Processing library: {}\n'.format(libID))
        KDE = _get_KDEs_for_libID(KDEs, KDE_type, libID, comm)

        # if needed: making a list of taxa that can incorporate 
        ## (unique to this library)
        ## TODO: abundance cutoff: 
        ### taxa must have abundance > threshold to incorporate
        ## TODO: abundance weighting with less incorp for less taxa
        try:
            incorp_list_len = len(taxa_incorp_list[libID])
        except KeyError:
            msg = 'Library "{}" not found'.format(libID)
            raise KeyError, msg
        if incorp_list_len == 0:
            taxa_incorp_list[libID] = _taxon_incorp_list(libID, config, KDE)
            
        # setting params for parallelized function
        pfunc = partial(_make_kde, 
                        config = config,
                        libID = libID,
                        stat_dir = tmpdirpath,
                        n = args['-n'],
                        taxa_incorp_list = taxa_incorp_list[libID],
                        isotope = args['--isotope'],
                        bw_method = args['--bw'])                    

        # parallel by taxon
        pool = ProcessingPool(nodes=int(args['--np']))
        if args['--debug']:
            tmp = map(pfunc, KDE.items())
        else:
            tmp = pool.map(pfunc, KDE.items())
        KDE = None
        # storing dict of KDEs 
        if args['-o'].lower() == 'none':            
            KDEs_iso[libID] = {taxon:kde for taxon,kde in tmp}
        else:
            KDEs_iso[libID] = Utils.write_lib_kde({taxon:kde for taxon,kde in tmp},
                                                  args['-o'], 
                                                  libID)
        tmp = None
    
    # combine stats
    _write_stats(args['--shift'], tmpdirpath)
    shutil.rmtree(tmpdirpath)

    # writing pickled BD-KDE with BD shift from isotope incorp
    if args['-o'].lower() == 'none':
        dill.dump(KDEs_iso, sys.stdout)    
    else:
        with open(args['-o'], 'wb') as outFH:
            dill.dump(KDEs_iso, outFH)


def _get_KDEs_for_libID(KDEs, KDE_type, libID, comm=None):
    """Parse out dict of KDE objects for just libID.
    Parsing depends on the KDE type
    """
    if KDE_type == 1:
        KDE = {t:k for t,k in KDEs}
    elif KDE_type == 2:
        KDE = KDEs
    elif KDE_type == 3:
        try: 
            KDE = KDEs[libID]
        except KeyError:
            ## kde library not found, duplicating KDE
            msg = 'WARNING: config library {} not found in KDEs.' + \
                  'Using a different KDE object\n'
            sys.stderr.write(msg.format(libID))
            KDE = KDEs[KDEs.keys()[0]]
    elif KDE_type == 4:
        try:
            KDE = Utils.load_kde(KDEs[libID])
        except KeyError:
            ## kde library not found, duplicating KDE
            msg = 'WARNING: config library {} not found in KDEs.' + \
                  'Using a different KDE object\n'
            sys.stderr.write(msg.format(libID))
            KDE = Utils.load_kde(KDEs[KDEs.keys()[0]])
            # adding comm info to KDEs 
        _add_comm_to_kde(KDE, comm)
    else:
        raise ValueError, 'KDE object type not recognized'
    return KDE


def _make_kde(x, libID, config, taxa_incorp_list, 
              isotope='13C', n=10000, bw_method=None, stat_dir=None): 
    """Making new KDE of BD value distribution which includes
    BD shift due to isotope incorporation. 

    Parameters
    ----------
    x : list
        [taxon_name, dict -- {kde:abundance}]
    libID : str
        library ID
    config : config object
    taxa_incorp_list : list
        names of taxa that can incorporate isotope
    isotope : str, optional
        isotope that is incorporated
    n : int
        number of Monte Carlo samples to use for estimating BD+isotope_BD
        distribution
    bw_method : str or function
        bandwidth scalar or function passed to scipy.stats.gaussian_kde().
    stat_dir : str
        directory path for writing BD shift stats. Nothing written if None
    Returns
    -------
    tuple -- (taxon_name, KDE*)
        *Note: KDE object may be None    
    """
    taxon_name,x = x
    try:
        bw_method = float(bw_method)
    except (TypeError, ValueError) as e:
        pass

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
        shift_stats = [libID, taxon_name] + [float(0)] * 6
        _write_tmp_stats(shift_stats, stat_dir)
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

    # BD shift stats
    BDs = kde.resample(n)[0]
    BDs_wIncorp = add_incorp(np.copy(BDs), mix_model, maxIsotopeBD)
    shift_stats = _calc_BD_shift_stats(libID, taxon_name, BDs, BDs_wIncorp)
    _write_tmp_stats(shift_stats, stat_dir)
    
    # making KDE of BD + BD_isotope_incorp
    kdeBD = stats.gaussian_kde(BDs_wIncorp, bw_method=bw_method)

    # return
    return (taxon_name, kdeBD)


def _calc_BD_shift_stats(libID, taxon_name, BDs, BDs_wIncorp):
    BD_shift = BDs_wIncorp - BDs
    BDs = None
    shift_stats = [libID, taxon_name, 
                   np.min(BD_shift), 
                   np.percentile(BD_shift, 25),
                   np.mean(BD_shift),
                   np.median(BD_shift), 
                   np.percentile(BD_shift, 75), 
                   np.max(BD_shift)]
    BD_shift = None
    return shift_stats


def __add_comm_to_kdes(taxon_name, kde, comm, libIDs):    
    d = {'kde':kde, 'abundances':{}}                                         
    for libID in libIDs:                                                        
        abund = comm.get_taxonAbund(taxon_name, libID=libID)
        try:
            d['abundances'][libID] = abund[0]
        except IndexError:
            msg = 'WARNING; no abundance data for: lib={}, taxon={}\n'
            sys.stderr.write(msg.format(libID, taxon_name))
            d['abundances'][libID] = 0
    return d


def _add_comm_to_kde(KDEs, comm):
    """Adding comm data for each taxon to each KDE.
    'abundances' will be an empty dict if comm is not provided.
    In-place edit of KDE_BD {taxon_name:{kde|abundances}}

    Parameters
    ----------
    KDE_BD : KDE object       
    comm : gradient community object
    """
    try:
        libIDs = comm.get_unique_libIDs()
    except AttributeError:
        libIDs = []

    for x,y in KDEs.items():        
        try: 
            d = {}
            for xx,yy in y.items():                
                d[xx] =  __add_comm_to_kdes(xx, yy, comm, libIDs)
            KDEs[x] = d
        except AttributeError:
            KDEs[x] = __add_comm_to_kdes(x, y, comm, libIDs)


def _write_tmp_stats(stats, dirpath):
    if dirpath is None:
        return 0
    outfn = tempfile.NamedTemporaryFile()
    outfn = os.path.split(outfn.name)[1]
    outfn = os.path.join(dirpath, outfn + '_stats.txt')
    stats = [str(x) for x in stats]
    with open(outfn, 'wb') as outfh:
        outfh.write('\t'.join(stats) + '\n')

    
def _write_stats(outfn, tmpdirpath):
    tmpfiles = glob.glob(os.path.join(tmpdirpath, '*_stats.txt'))
    if len(tmpfiles) == 0:
        return 0
    header = ['library', 'taxon', 'min', 'q25', 'mean', 'median', 
              'q75', 'max']
    df = [] 
    for F in tmpfiles:
        with open(F, 'rb') as infh:
            for line in infh:
                line = line.rstrip().split('\t')
                df.append(line)

    df = pd.DataFrame(df, columns=header)
    #df = df.sort(['library', 'taxon'])
    df = df.sort_values(by=['library','taxon'])
    df.to_csv(outfn, sep='\t', index=False)
    sys.stderr.write('File written: {}\n'.format(outfn))


def _taxon_incorp_list(libID, config, KDE_BD):
    """Make a list of taxa that incorporated isotope.
    """
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

    Parameters
    ----------
    taxon_name : str
        taxon name string
    libID : str
        library ID string
    config : config object

    Returns
    -------
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
    """Select the intra-population parameter value
    based on the inter-population distribution function.
    Values are % isotope incorporation, so acceptable
    range is 0-100 (will try 'maxtries' times to select value in range).

    Parameters
    ----------
    interPopDist : dict
        {'interPopDist':{'function':interPopdist_function}}
    taxon_name : str
         name of taxon
    maxtries : int
         number of tries to get a parameter values >0

    Returns
    -------
    float : intra-pop param value
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
    In-place edit of params.

    Parameters
    ----------
    params : dict
        {param_ID:param_value}
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
    """Setting the theoretical max BD shift of an isotope 
    (if 100% incorporation).

    Parameters
    ----------
    isotope : str
        name of isotope

    Returns
    -------
    float : max BD value
    """
    psblIsotopes = {'13C' : 0.036,
                    '15N' : 0.016}
    try:
        return psblIsotopes[isotope.upper()]
    except KeyError:
        raise KeyError('Isotope "{}" not supported.'.format(isotope))


def _load_taxa_incorp_list(inFile, config):
    """Loading list of taxa that incorporate isotope.
    Parameters
    ----------
    inFile : str
        File name of taxon list
    config : config object

    Returns
    -------
    {library:[taxon1, ...]}
    """
    taxa = {}
    with open(inFile, 'rb') as inFH:
        for line in inFH:
            line = line.rstrip().split('\t')
            # if 1 column, using config-defined libraries
            if len(line) == 1:
                line = [[x,line[0]] for x in config.keys()]
            else:
                line = [line]
            
            for x in line:
                try:
                    taxa[x[0]].append(x[1])
                except KeyError:
                    taxa[x[0]] = [x[1]]

    return taxa

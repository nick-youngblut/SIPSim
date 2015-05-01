# import
## batteries
import os, sys
import re
import functools
import types
import inspect
import logging
from collections import defaultdict
from random import shuffle
from StringIO import StringIO
## 3rd party
import configobj
from configobj import ConfigObj, flatten_errors
from validate import Validator
import numpy as np
import pandas as pd
import mixture
## application
import Utils
import TraitEvo
from TraitEvo import BrownianMotion
import SIPSimCython


# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def get_configspec(strIO=True):
    """Return configspec set for instance.
    Args:
    strIO -- return configspec as a StringIO instance
    """    
    configspec = """
    [__many__]
       max_perc_taxa_incorp = float(0,100, default=100)

       [[__many__]]
          distribution = option('normal','uniform', 'BM', default='normal')
          weight = float(0,1, default=1)
    
            [[[__many__]]]
                [[[[__many__]]]]
                distribution = option('normal','uniform', 'BM', default='normal')
                weight = float(0,1, default=None)
                mu = float(0,100, default=None)
                sigma = float(0,100, default=None)
                start = float(0,100, default=None)
                end = float(0,100, default=None)  
                ratio = float(0,1, default=None)
                minVal = float(0,100, default=None)
                maxVal = float(0,100, default=None)
    """
    if strIO == True:
        return StringIO(configspec)
    else:
        return configspec
    

def make_incorp_model(taxon_name, libID, config):
    """Setting isotope incorporation based on the interPopDist
    function for each intraPop parameter.
    Args:
    taxon_name -- taxon name string
    libID -- library ID string
    config -- config object
    Return:
    mixture model object (mixture class)
    """
    psblDists = {'normal' : mixture.NormalDistribution,
                 'uniform' : mixture.UniformDistribution}

    # creating individual distribution functions
    libSect = config.get_libSection(libID[0])
    intraPopDist_IDs = []
    intraPopDist_funcs = []
    weights = []
    for (intraPopDistID,intraPopDist) in config.iter_configSections(libSect):    
        # name of standard distribution
        distID = intraPopDist['distribution']
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
    weights = _fill_in_weights(weights)
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
    Return:
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
    Return:
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


def _fill_in_weights(weights, n=None, total=1.0):
    """Filling in any missing weights (None) 
    so that the sum of all weights equals 'total'.
    Both None and zero values will be set to 0,
    while all others scaled to sum to 'total'.
    Args:
    weights -- iterable of weight values
    n -- number of total weights needed.
         if None; len(weights) used.
    total -- sum of all returned weight values
    Return:
    list of floats
    """
    total = float(total)
        
    # adding Nones if needed
    if n is None:
        n = len(weights)
    elif n > len(weights):
        i = n - len(weights)
        weights.append([None for x in xrange(i)])
    elif n < len(weights):
        msg = 'WARNING: "n" is < len(weights)\n'
        sys.stderr.write(msg)
    else:
        pass

    # summing all real values
    w_sum = sum([x for x in weights if x is not None])


    # setting weights to sum to total
    ## scaling factor
    if w_sum == 0:
        return [total / n for x in xrange(n)]
    else:
        a = total / w_sum

    ## scaling weights
    if a == 1:
        # convert Nones to zeros
        weights= [x if x else 0 for x in weights]
    elif a > 1:
        # convert Nones to zeros
        weights = [x if x else 0 for x in weights]
        # rescaling to total 
        weights = [x * a for x in weights]
    elif 1 < 1:
        # rescaling all real values to sum to total
        # all None values set to 0
        ## rescale
        weights = [x * a if x else x for x in weights]
        ## convert Nones to zeros
        weights = [x if x else 0 for x in weights]    
    else:
        raise RuntimeError, 'Logic Error'

    assert sum(weights) == total, 'Scaled weights != total'
    return weights


def isotopeMaxBD(isotope):
    """Setting the theoretical max BD shift of an isotope (if 100% incorporation).
    Args:
    isotope -- str; name of isotope
    Return:
    float
    """
    psblIsotopes = {'13C' : 0.036,
                    '15N' : 0.016}
    try:
        return psblIsotopes[isotope.upper()]
    except KeyError:
        raise KeyError('Isotope "{}" not supported.'.format(isotope))



        
def populationDistributions_OLD(config, comm):
    """Setting distribution parameters for intra-population isotope incorporation
    Args:
    config -- IsoIncorpConfig instance
    comm -- CommTable instance
    """    
    # list of taxon names used for all library iterations
    taxon_names = comm.get_unique_taxon_names()
    ## TODO: introduced biased shuffling based on taxon abundance/rank?
    shuffle(taxon_names)
    
    # DataFrame for output
    outTbl = pd.DataFrame(columns=['library', 'taxon_name', 
                                   'distribution_index', 'distribution', 
                                   'weight', 'param', 'param_value'])


    # iterating through all libraries & taxon in each library
    configLibIDs = {k:v for k,v in zip(config.get_configLibs(trunc=True), \
                                       config.get_configLibs())}
    for libID in comm.iter_libraries():        
        # cross check libID between comm and config files
        libID = str(libID)
        if not libID in configLibIDs:
            msg = 'library "{}" in comm file but not in config. Skipping\n'
            logging.warning(msg.format(libID))
            continue
        else:
            libID = [x for x in (configLibIDs[libID], libID)]  # (full,truncated)
            libSect = config.get_libSection(libID[0])

        # max percent of taxa that can have isotope incorporation (set in config)
        max_perc_inc = config.get_max_perc_taxa_incorp(libID[0])

        # iterating by taxon
        ntaxa = len(taxon_names)
        for taxon_idx in xrange(ntaxa):
            taxon_name = taxon_names[taxon_idx]

            # check if taxon in library
            if not comm.taxonInLib(taxon_name, libID[0]):
                continue                

            # writing out the intra-pop dist(s) for the taxon in the library
            ## only taxa in the user-defined % should have any incorporation
            perc_taxa_processed = float(taxon_idx + 1) / ntaxa * 100
            if(perc_taxa_processed > max_perc_inc):            
                _set_intraPopDistZero(outTbl, libID[1], taxon_name)
            else:
                _set_intraPopDist(outTbl, config, libID, taxon_name)
                        
    return outTbl

                
def _set_intraPopDist_OLD(outTbl, config, libID, taxon_name):
    """Setting isotope incorporation based on the interPopDist
    function for each intraPop parameter.
    Args:
    outTbl -- output dataframe
    config -- config 
    libID -- library ID [full_version,truc_version]
    taxon_name -- taxon name
    """

    libSect = config.get_libSection(libID[0])
    dist_cnt = 0
    for (intraPopDistID,intraPopDist) in config.iter_configSections(libSect):    
        dist_cnt += 1
        # getting distribution
        distribution = intraPopDist['distribution']
        # getting weight
        weight = float(intraPopDist.get('weight', 1.0))
            
        # setting intra-pop dist param values
        intraPopParams = dict()
        for (paramID,param) in config.iter_configSections(intraPopDist):
            # getting inter-pop dist function
            func = [v['function'] for k,v in config.iter_configSections(param)][0]

            # sampling from function to get parameter for intra-pop distribution
            maxtries = 1000
            tries = 0
            while True:
                tries += 1
                try:
                    paramVal = func.sample(taxon_name)
                except TypeError:
                    paramVal = func.sample()
                try:
                    paramVal = paramVal[0]
                except TypeError:
                    pass
                # values must be in range: 0-100
                if paramVal >= 0 and paramVal <= 100:
                    intraPopParams[paramID] = paramVal
                    break
                # exceeding maxtries?
                if tries >= maxtries:
                    er = 'Taxon: {}'.format(taxon_name)
                    msg = 'Exceeded maxtries to get parameter in range: 0-100'  
                    sys.exit(': '.join(err, msg))

        # checking paramters
        _check_intraPopParams(intraPopParams)
        
        # writing values
        for paramID, paramVal in intraPopParams.items():
            outTbl.loc[outTbl.shape[0]+1] = [libID[1], taxon_name, 
                                             str(dist_cnt), distribution,
                                             str(weight), paramID, str(paramVal)]


def _check_intraPopParams_OLD(intraPopParams):
    """Formatting and verifying intra-population parameters.
    Args:
    intraPopParams -- dict {param : paramValue}
    Return:
    None -- modifies dict in place
    """
    # mu values rounded if < 1e-10
    try:
        if intraPopParams['mu'] < 1e-10:
            intraPopParams['mu'] = round(intraPopParams['mu'], 0)
    except KeyError:
        pass
            
    # if start-end params differ by < 1e-8, make same value
    try:
        if abs(intraPopParams['start'] - intraPopParams['end']) < 1e-8:
            intraPopParams['start'] = round(intraPopParams['start'],0)
            intraPopParams['end'] = round(intraPopParams['end'],0)
    except KeyError:
        pass

    # flip start-end if needed
    try:
        if intraPopParams['start'] > intraPopParams['end']:
            intraPopParams['end'],intraPopParams['start'] = (intraPopParams['start'],
                                                                 intraPopParams['end'])
    except KeyError:
        pass   
                                                         
    return None

                                                    
def _set_intraPopDistZero_OLD(outTbl, libID, taxon_name):
    """Setting isotope incorporation essentially to 0
    (distribution=uniform, start=0, end=0)
    Args:
    outTbl -- output pandas dataframe
    libID -- library ID
    taxon_name -- taxon name
    """
    # saving row of values
    l = [libID, taxon_name, '1', 'uniform', '1.0']
    outTbl.loc[outTbl.shape[0]+1] = l + ['start', '0.0']
    outTbl.loc[outTbl.shape[0]+1] = l + ['end', '0.0'] 


        
class Config(ConfigObj):
    """Subclassing ConfigObj for isotope incorporation file."""

    @classmethod
    def load_config(cls, config_file, phylo=None):
        configspec = get_configspec()

        # loading config file
        return cls(config_file, configspec=configspec, phylo=phylo)


    def __init__(self, *args, **kwargs):

        # phylo
        self.phylo = kwargs.pop('phylo', None)
        self.tree = TraitEvo.read_tree(self.phylo)

        # config init
        ConfigObj.__init__(self, *args, **kwargs)
        
        # validate
        self._validate_config()
        
        # check param assertions
        self._param_asserts()

        # setting inter-population distribution functinos
        self._set_interPopDistFuncs()
        self._set_interPopDistMM()
                

    def _validate_config(self):
        """Run configobj validator on config.
        Args:
        config -- ConfigObj instance
        """
        vtor = Validator()
        res = self.validate(vtor, preserve_errors=True)
        if res != True:
            print 'ERROR: config file contains errors:'
            
            for entry in flatten_errors(self, res):
                # each entry is a tuple
                section_list, key, error = entry
                if key is not None:
                    section_list.append(key)
                else:
                    section_list.append('[missing section]')
                    section_string = ', '.join(section_list)
                if error == False:
                    error = 'Missing value or section.'
                print section_string, ' = ', error
            sys.exit(1)

            
    def _param_asserts(self):
        """Checking that certain parameters comply to assertions.
        """    
        for (libID,sect1) in self.iteritems():
            for (intraPopDist,sect2) in self.iter_configSections(sect1):
                # intra-pop level
                self._check_start_end(sect2)
                for (intraPopDistParam,sect3) in self.iter_configSections(sect2):
                    for (interPopDistParam,sect4) in self.iter_configSections(sect3):
                        # inter-pop level
                        self._check_start_end(sect4)


    def _check_start_end(self, sect):
        """Checking start & end parameters.
        Args:
        sect -- config section class
        """
        keyParams = {k.lower():v for k,v  in self.iter_sectionKeywords(sect)}
                       
        # assert params: end  > start
        if ('start' in keyParams) & ('end' in keyParams):
            try:
                startVal = float(keyParams['start'])
                endVal = float(keyParams['end'])
            except TypeError:
                return None
                
            if startVal > endVal:
                sect['start'] = endVal
                sect['end'] = startVal                
            elif startVal == endVal:                
                if endVal >= 100:
                    sect['start'] = float(sect['start']) - 1e-10
                else:
                    sect['end'] = float(sect['end']) + 1e-10
            else:
                pass

        return None
                                            

    def _set_interPopDistFuncs(self):
        """Setting the inter-population distribution random sampling functions
        based on user-selected distributions (pymix distributions).
        A mixture models will be used for the inter-pop distribution
        if multiple distributions are provided
        (e.g., [[[[interPopDist 1]]]] & [[[[interPopDist 2]]]]).
        'weight' parameters can be set for the individual distributions;
        otherwise, each individual distribution will be weighted equally. 

        """
        psblDists = {'normal' : mixture.NormalDistribution,
                     'uniform' : mixture.UniformDistribution,
                     'BM' : BrownianMotion}

        for (libID,sect1) in self.iteritems():
            for (intraPopDist,sect2) in self.iter_configSections(sect1):
                for (intraPopDistParam,sect3) in self.iter_configSections(sect2):        
                    for (interPopDist,sect4) in self.iter_configSections(sect3):

                        dist = sect4['distribution']                 
                        otherParams = {k:v for k,v in sect4.items() \
                                       if k not in ['distribution','weight'] \
                                       and v is not None}

                        try: 
                            sect4['function'] = psblDists[dist](tree=self.tree, **otherParams)
                        except TypeError:
                            try:
                                sect4['function'] = psblDists[dist](**otherParams)            
                            except KeyError:
                                msg = 'distribution "{}" not supported'
                                raise KeyError(msg.format(dist))


    def _set_interPopDistMM(self):
        """Setting the inter-population distribution random sampling functions
        based on user-selected distributions (pymix distributions).
        A mixture model will be used for the inter-pop distribution
        if multiple distributions are provided
        (e.g., [[[[interPopDist 1]]]] & [[[[interPopDist 2]]]]).
        'weight' parameters can be set for the individual distributions;
        otherwise, each individual distribution will be weighted equally. 
        """
        for (libID,sect1) in self.iteritems():
            for (intraPopDist,sect2) in self.iter_configSections(sect1):
                for (intraPopDistParam,sect3) in self.iter_configSections(sect2):                
                    
                    # setting mixture models 
                    allInterPopDists = [[k,v] for k,v in self.iter_configSections(sect3)]
                    nInterPopDists = len(allInterPopDists)

                    allFuncs = [v['function'] for k,v in allInterPopDists]
                    allWeights = [v['weight'] for k,v in allInterPopDists \
                                  if 'weight' in v]

                    BMfuncs = [x for x in allFuncs if isinstance(x, TraitEvo.BrownianMotion)]
                    
                    # no mixture model if BM is a selected distribution
                    if len(BMfuncs) >= 1:
                        sect3['interPopDist'] = {'function' : BMfuncs[0]}
                    else:

                        # else create mixture model from >=1 standard distribution
                        ## fill in missing weights (total = 1)
                        allWeights = _fill_in_weights(allWeights, n = nInterPopDists) 
                        assert sum(allWeights) == 1, 'interpop weights don\'t sum to 1'
                        
                        # changing interPopDist to mixture model
                        for i in sect3: 
                            sect3.popitem()                        
                        sect3['interPopDist'] = {'function' :
                                                 mixture.MixtureModel(nInterPopDists,
                                                                      allWeights,
                                                                      allFuncs)}

            

    def get_max_perc_taxa_incorp(self, libID):
        """Getting the user-defined max_perc_taxa_incorp
        for a specified library.
        Args:
        libID -- library ID
        Return:
        float
        """
        try:
            return self[libID]['max_perc_taxa_incorp']
        except KeyError:
            msg = 'Cannot find max_perc_taxa_incorp for library: {}'
            raise KeyError, msg.format(libID)

                    
    def get_libSection(self, libID):
        """Getting sub-section of user-defined library from config.
        Args:
        libID -- library ID
        """
        return self[libID]

        
    def get_configLibs(self, trunc=False):
        """Get the library IDs from the config
        """
        if trunc==True:
            return [re.sub('\D+', '', x) for x in self.keys()]
        else:
            return self.keys()

            
    def iter_configSections(self, sect=None):
        """Iter through the sub-sections of the provided section.
        If not section provided, starting at the top-level section.
        """
        if sect is None:
            sect = self
        for (sectID,val) in sect.iteritems():
            if isinstance(val, configobj.Section):
                yield (sectID,val)
                            
                
    def iter_sectionKeywords(self, sect=None):
        """Iter through just the keywords of the provided section.
        If not section provided, starting at the top-level section.        
        """
        if sect is None:
            sect = self
        for (sectID,val) in sect.iteritems():
            if not isinstance(val, configobj.Section):
                yield (sectID,val)

                

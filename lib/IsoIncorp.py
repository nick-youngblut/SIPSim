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


# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def get_configspec(strIO=True):
    """Return configspec set for instance.
    Args:
    strIO -- return configspec as a StringIO instance
    """    
    configspec = """
    [__many__]
       max_percent_incorp = float(0,100, default=100)

       [[__many__]]
          distribution = option('normal','uniform', 'BM', default='normal')
          weight = float(0,1, default=None)
    
            [[[__many__]]]
                [[[[__many__]]]]
                distribution = option('normal','uniform', default='normal')
                weight = float(0,1, default=None)
                mu = float(0,100, default=None)
                sigma = float(0,100, default=None)
                start = float(0,100, default=None)
                end = float(0,100, default=None)    
    """
    if strIO == True:
        return StringIO(configspec)
    else:
        return configspec

        
def populationDistributions(config, comm):
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
    outTbl = pd.DataFrame(columns=['library', 'taxon_name', 'distribution_index',
                                   'distribution', 'weight', 'param', 'param_value'])


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

        # max percent taxa that can have isotope incorporation
        max_perc_inc = config.get_max_percent_incorp(libID[0])

        # iterating by taxon
        ntaxa = len(taxon_names)
        for taxon_idx in xrange(ntaxa):
            taxon_name = taxon_names[taxon_idx]

            # check if taxon in library
            if not comm.taxonInLib(taxon_name, libID[0]):
                continue                

            # writing out the intra-population distribution(s) for the taxon in the library
            ## only taxa in the user-defined percentage should have any incorporation
            perc_taxa_processed = float(taxon_idx + 1) / ntaxa * 100
            if(perc_taxa_processed > max_perc_inc):            
                _set_intraPopDistZero(outTbl, libID[1], taxon_name)
            else:
                _set_intraPopDist(outTbl, config, libID, taxon_name)
                        
    return outTbl

                
def _set_intraPopDist(outTbl, config, libID, taxon_name):
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
                    msg = 'Exceeded maxtries to get parameter in range: 0-100'
                    sys.exit(msg)

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
        
        # writing values
        for paramID, paramVal in intraPopParams.items():
            outTbl.loc[outTbl.shape[0]+1] = [libID[1], taxon_name, 
                                             str(dist_cnt), distribution,
                                             str(weight), paramID, str(paramVal)]

                                                    
def _set_intraPopDistZero(outTbl, libID, taxon_name):
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

    def __init__(self, *args, **kwargs):
        # phylo
        self.phylo = kwargs.pop('phylo', None)

        # config init
        ConfigObj.__init__(self, *args, **kwargs)
        
        # validate
        self._validate_config()
        
        # check param assertions
        self._param_asserts()

        # setting inter-population distribution functinos
        self._set_interPopDistFuncs()
        self._set_interPopDistMM()
        #self._set_intraPopDistFuncs()


    @classmethod
    def load_config(cls, config_file, phylo=None):
        configspec = get_configspec()

        # loading config file
        return cls(config_file, configspec=configspec)
    

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
                
            if startVal >= endVal:
                if endVal >= 100:
                    sect['start'] = float(sect['start']) - 1e-10
                else:
                    sect['end'] = float(sect['end']) + 1e-10
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
        # TODO: add BM distribution
        psblDists = {'normal' : mixture.NormalDistribution,
                     'uniform' : mixture.UniformDistribution}

        for (libID,sect1) in self.iteritems():
            for (intraPopDist,sect2) in self.iter_configSections(sect1):
                for (intraPopDistParam,sect3) in self.iter_configSections(sect2):        
                    for (interPopDist,sect4) in self.iter_configSections(sect3):

                        dist = sect4['distribution']                 
                        otherParams = {k:v for k,v in sect4.items() \
                                       if k not in ['distribution','weight'] \
                                       and v is not None}

                        try:
                            sect4['function'] = psblDists[dist](**otherParams)            
                        except KeyError:
                            msg = 'distribution "{}" not supported'
                            raise KeyError(msg.format(dist))


    def _set_interPopDistMM(self):
        """Setting the inter-population distribution random sampling functions
        based on user-selected distributions (pymix distributions).
        A mixture models will be used for the inter-pop distribution
        if multiple distributions are provided
        (e.g., [[[[interPopDist 1]]]] & [[[[interPopDist 2]]]]).
        'weight' parameters can be set for the individual distributions;
        otherwise, each individual distribution will be weighted equally. 
        """
        for (libID,sect1) in self.iteritems():
            for (intraPopDist,sect2) in self.iter_configSections(sect1):
                for (intraPopDistParam,sect3) in self.iter_configSections(sect2):        

                    # setting mixture models if needed
                    allInterPopDists = [[k,v] for k,v in self.iter_configSections(sect3)]
                    nInterPopDists = len(allInterPopDists)
                    if nInterPopDists > 1:
                        allFuncs = [v['function'] for x in allInterPopDists]
                        allWeights = [v['weight'] for k,v in allInterPopDists \
                                      if 'weight' in v]
                        # fill in missing weights (total = 1)
                        allWeights = self._fill_in_weights(allWeights, n = nInterPopDists) 
                        assert sum(allWeights) == 1, 'interpop weights don\'t sum to 1'
                            
                        # changing interPopDist to mixture model
                        for i in sect3: sect3.popitem()                        
                        sect3['interPopDist 1'] = {'function' :
                                                   mixture.MixtureModel(nInterPopDists,
                                                                        allWeights,
                                                                        allFuncs)}
                                        

    def _fill_in_weights(self, allWeights, n=None, total=1):
        """Filling in any missing weights so that all weights total to 'total'
        Args:
        allWeights -- iterable of weight values
        n -- number of weight values returned (if None, uses len(allWeights)
        total -- sum of all returned weight values
        """
        if n is None:
            n = len(allWeigths)

        # convert NoneTypes
        allWeightsWithVals = [x for x in allWeights if x is not None]
        cur_sum = sum(allWeightsWithVals)
        
        if cur_sum > total:
            raise ValueError, "Sum of weights > 1"
        elif len(allWeightsWithVals) > n:
            raise ValueError, "len(allWeights) > n"
        elif cur_sum == total and len(allWeightsWithVals) == n:
            return allWeightsWithVals
        else:
            remainder = float(total - cur_sum)
            nNeeded = n - len(allWeightsWithVals)
            newWeights = [remainder / nNeeded] * nNeeded
            for i,val in enumerate(allWeights):
                if val is None:
                    allWeights[i] = newWeights.pop()
            return allWeights + newWeights
            

    def get_max_percent_incorp(self, libID):
        """Getting the user-defined max_percent_incorp
        for a specified library.
        Args:
        libID -- library ID
        """
        try:
            return self[libID]['max_percent_incorp']
        except KeyError:
            msg = 'Cannot find max_percent_incorp for library: {}'
            raise KeyError, msg.format(libID)

                    
    def get_libSection(self, libID):
        """Getting sub-section of user-defined library from config.
        Args:
        libID -- library ID
        """
        return self[libID]

        
    def get_configLibs(self, trunc=False):
        """Get the library IDs from the config"""
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

                


# class same_value(object):

#     def __init__(self, start, end, **kwargs):
#         self.start = start
#         self.end = end
        
#     def p(self):
#         return self.start

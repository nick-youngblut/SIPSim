# import
## batteries
import os, sys
import re
import functools
import types
from collections import defaultdict
from random import shuffle
import inspect
from StringIO import StringIO
import logging
## 3rd party
import configobj
from configobj import ConfigObj, flatten_errors
from validate import Validator
import numpy as np
#import pymix.mixture as mixture
import mixture
import pandas as pd
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

        
def populationDistributions(config, comm, percTaxa=100):
    """Setting distribution parameters for intra-population isotope incorporation
    Args:
    config -- IsoIncorpConfig instance
    comm -- CommTable instance
    percTaxa -- percent of taxa to have any isotope incorporation
    """
    # assert args
    percTaxa = float(percTaxa)
    
    # list of taxon names used for all library iterations
    taxon_names = comm.get_unique_taxon_names()
    ## TODO: introduced biased shuffling based on taxon abundance/rank?
    shuffle(taxon_names)
    
    # DataFrame for output
    outTbl = pd.DataFrame(columns=['library', 'taxon_name', 'distribution_index',
                                   'distribution', 'weight', 'param', 'param_value'])


    # iterating through all libraries & taxon in each library
    configLibIDs = {k:v for k,v in zip(config.get_configLibs(trunc=True), config.get_configLibs())}
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

        # iterating by taxon
        ntaxa = len(taxon_names)
        for taxon_idx in xrange(ntaxa):
            taxon_name = taxon_names[taxon_idx]

            # check if taxon in library
            if not comm.taxonInLib(taxon_name, libID[0]):
                continue
                
            # writing out the intra-population distribution(s) for the taxon in the library
            if((float(taxon_idx)+1) / ntaxa * 100 > percTaxa):            
                set_intraPopDistZero(outTbl, libID[1], taxon_name)
            else:
                set_intraPopDist(outTbl, config, libID, taxon_name)
                        
    return outTbl

                
def set_intraPopDist(outTbl, config, libID, taxon_name):
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
        try:
            weight = intraPopDist['weight']
        except KeyError:
            weight = ''
            
        # setting intra-pop param values
        for (paramID,param) in config.iter_configSections(intraPopDist):
            # getting inter-pop dist function
            func = [v['function'] for k,v in config.iter_configSections(param)][0]
            # sampling from function to get parameter for intra-pop distribution
            while True:
                paramVal = func.sample()
                try:
                    paramVal = paramVal[0]
                except TypeError:
                    pass
                # values must be in range: 0-100
                if paramVal >= 0 and paramVal <= 100:
                    break

            # setting weight as float
            if weight is None:
                weight = 1.0
            else:
                weight = float(weight)
                    
            # saving row of values
            outTbl.loc[outTbl.shape[0]+1] = [libID[1], taxon_name, str(dist_cnt), distribution,
                                             str(weight), paramID, str(paramVal)]

                                                    
def set_intraPopDistZero(outTbl, libID, taxon_name):
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
            sys.exit()

            
    def _param_asserts(self):
        """Checking that certain parameters comply to assertions.
        """
        # TODO: finish
        
        def end_gt_start(x):
            pass
    
        for (libID,val1) in self.iteritems():
            for (intraPopDist,val2) in self.iter_configSections(val1):
                for (intraPopDistParam,val3) in self.iter_configSections(val2):
                    for (interPopDistParam,val4) in self.iter_configSections(val3):

                        keyParams = {k.lower():v for k,v  in self.iter_sectionKeywords(val4)}
                        
                        # assert params: end  > start
                        if ('start' in keyParams) & ('end' in keyParams):
                            try:
                                startVal = float(keyParams['start'])
                                endVal = float(keyParams['end'])
                            except TypeError:
                                continue
                                
                            if startVal >= endVal:
                                if endVal >= 100:
                                    val4['start'] = float(val4['start']) #- 1e-10
                                else:
                                    val4['end'] = float(val4['end']) #+ 1e-10
                                    
                                            

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
                     'uniform' : mixture.UniformDistribution,
                     'same_value' : lambda start,end: start
                 }

        for (libID,val1) in self.iteritems():
            for (intraPopDist,val2) in self.iter_configSections(val1):
                for (intraPopDistParam,val3) in self.iter_configSections(val2):
                    # setting standard distributions
                    for (interPopDist,val4) in self.iter_configSections(val3):
                        d = val4['distribution']
                        otherParams = {k:v for k,v in val4.items() if k != 'distribution' and v is not None}

                        try: 
                            startParam = otherParams['start']
                            endParam = otherParams['end']
                        except KeyError:
                            pass
                        else:
                            # assert start <= end
                            if startParam > endParam:
                                otherParams['start'] = endParam
                                otherParams['end'] = startParam
                            elif startParam == endParam:  
                                if startParam >= 100:
                                    otherParams['start'] -= 1e-5
                                else:
                                    otherParams['end'] += 1e-5

                        
                        try:
                            val4['function'] = psblDists[d](**otherParams)
                        except KeyError:
                            raise KeyError('distribution "{}" not supported'.format(d))

                            
                    # setting mixture models if needed
                    allInterPopDists = [[k,v] for k,v in self.iter_configSections(val3)]
                    nInterPopDists = len(allInterPopDists)
                    if nInterPopDists > 1:
                        allFuncs = [v['function'] for k,v in allInterPopDists]
                        allWeights = [v['weight'] for k,v in allInterPopDists if 'weight' in k]
                        if len(allWeights) < nInterPopDists:
                            allWeights = [1.0 / nInterPopDists for i in xrange(nInterPopDists)]

                        # changing interPopDist to mixture model
                        for i in val3: val3.popitem()                        
                        val3['interPopDist 1'] = {'function' :
                                                  mixture.MixtureModel(nInterPopDists,
                                                                       allWeights,
                                                                       allFuncs)}

                    
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

                


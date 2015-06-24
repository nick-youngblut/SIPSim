# import
## batteries
import os, sys
import re
import functools
import types
import logging
from StringIO import StringIO
## 3rd party
import configobj
from configobj import ConfigObj, flatten_errors
from validate import Validator
import numpy as np
import pandas as pd
import mixture
## application
import TraitEvo
from TraitEvo import BrownianMotion


# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)



def get_configspec(strIO=True):
    """Return configspec set for instance.
    Args:
    strIO -- return configspec as a StringIO instance
    Returns:
    configspec object -- defines default parameters for config
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


def get_basicConfig(strIO=True):
    basicConfig = """
    [1]
      # baseline: no incorporation
    
      [[intraPopDist 1]]
      distribution = uniform
    
        [[[start]]]
    
          [[[[interPopDist 1]]]]
            distribution = uniform
            start = 0
            end = 0
    
        [[[end]]]
    
          [[[[interPopDist 1]]]]
            distribution = uniform
            start = 0
            end = 0
    
    [2]
      # 'treatment' community: possible incorporation
    
      max_perc_taxa_incorp = 50
    
      [[intraPopDist 1]]
      distribution = uniform
    
        [[[start]]]
    
          [[[[interPopDist 1]]]]
            distribution = uniform
            start = 100
            end = 100
    
        [[[end]]]
    
          [[[[interPopDist 1]]]]
            distribution = uniform
            start = 100
            end = 100
    """
    if strIO == True:
        return StringIO(basicConfig)
    else:
        return basicConfig


class ExampleConfig(ConfigObj):

    def __init__(self, *args, **kwargs):
        # config init
        ConfigObj.__init__(self, *args, **kwargs)
        
    def check_treatment_lib(self):
        try:
            self['2']
        except KeyError:
            self['2'] = {}    

    def check_treatment_intraPop(self):
        try:
            self['2']['intraPopDist 1'] 
        except KeyError:
            self['2']['intraPopDist 1'] = {}
                
    def set_percTaxa(self, percTaxa):
        self.check_treatment_lib()
        self['2']['max_perc_taxa_incorp'] = percTaxa

    def set_percIncorpUnif(self, percIncorpUnif):
        self.check_treatment_lib()
        self.check_treatment_intraPop()
        if self.treatment_intraPopDistDist != 'uniform':
            self.treatment_intraPopDistDist = 'uniform'
        self['2']['intraPopDist 1']['start'] = {
            'interPopDist 1' : {
                'distribution' : 'uniform',
                'start' : percIncorpUnif,
                'end' : percIncorpUnif
            }
        }        
        self['2']['intraPopDist 1']['end'] = {
            'interPopDist 1' : {
                'distribution' : 'uniform',
                'start' : percIncorpUnif,
                'end' : percIncorpUnif
            }
        }
        
    def set_percIncorpNorm(self, percIncorpMean, percIncorpSD):
        self.check_treatment_lib()
        self.check_treatment_intraPop()
        if self.treatment_intraPopDistDist != 'normal':
            self.treatment_intraPopDistDist = 'normal'
        for x in ['start', 'end']:
            try:
                self['2']['intraPopDist 1'].pop(x)
            except KeyError:
                pass

        self['2']['intraPopDist 1']['mu'] = {
            'interPopDist 1' : {
                'distribution' : 'uniform',
                'start' : percIncorpMean,
                'end' : percIncorpMean
            }
        }        
        self['2']['intraPopDist 1']['sigma'] = {
            'interPopDist 1' : {
                'distribution' : 'uniform',
                'start' : percIncorpSD,
                'end' : percIncorpSD
            }
        }
        

    @property
    def treatment_intraPopDistDist(self):
        try:
            return self['2']['intraPopDist 1']['distribution']
        except KeyError:
            msg = 'Cannot find "interPopDist 1" => "distribution"'
            raise KeyError, msg
        
    @treatment_intraPopDistDist.setter
    def treatment_intraPopDistDist(self, value):
        self['2']['intraPopDist 1']['distribution'] = value



class Config(ConfigObj):
    """Subclassing ConfigObj for isotope incorporation file.
    """
    @staticmethod
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
        Returns:
        list of floats -- [weight1, ...]
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



    @classmethod
    def load_config(cls, config_file, phylo=None):
        """Loading config object with validation via the pre-set configspec
        """
        configspec = get_configspec()

        cfg = cls(config_file, configspec=configspec, phylo=phylo)

        # checking keys
        msg = 'Library IDs must be a continuous set of integers starting from 1'
        def is_int(x):
            try: 
                return int(x)
            except ValueError:
                raise ValueError, msg
            return x
    
        libID_int = [is_int(x) for x in sorted(cfg.keys())]
        assert min(libID_int) == 1, msg
        assert sum(libID_int) == sum(range(len(libID_int)+1)), msg
        
        # return
        return cfg


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
        Returns:
        in-place edit
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
        Returns:
        None
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
        Returns:
        None
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
        Returns:
        None
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
                        allWeights = Config._fill_in_weights(allWeights, n = nInterPopDists) 
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
        Returns:
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
        try:
            return self[libID]
        except KeyError:
            msg = 'Cannot find library: "{}"'
            raise KeyError, msg.format(libID)
        
        
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

                

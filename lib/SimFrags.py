"""Classes and subclasses for genome fragment simulation"""

# import
## batteries
import os
import sys
import random
import functools

## 3rd party
import numpy as np
import pymc.distributions as pymcDist
import skew_normal



class SimFrags(object):
    """Class for genome fragment simulation"""

    def __init__(self, fld, flr, rtl=None):
        """
        Args:
        fld -- fragment length distribution (type,moments...)
        flr -- fragment length range (min,max)
        rtl -- read template length distribution (type,moments..)
        """
        
        self.flr = [float('inf') if x.lower() == 'none' or x.lower == 'inf' or x is None else int(x) for x in flr]

        # setting fld distribution function
        self.fld = self._set_user_dist(fld)
        assert self.fld is not None, 'fragment length distribution (fld) should not be None'
        
        # setting rtl distribution function
        if rtl is not None:
            self.rtl = self._set_user_dist(rtl)

        
    def simFrag(self, genome):
        """Simulation of a single fragment.
        Args:
        genome -- genome-like object
        Return:
        list -- [scaffold,start,end,frag_sequence]
        """

        tryCnt = 0
        while 1:
            # get read template position
            try:
                readTempPos = self._get_randPos_amp(genome)
            except AttributeError:
                readTempPos = self._get_randPos_shotgun(genome)
            assert len(readTempPos) == 3, 'readTempPos should be: "scaf,start,end"'
            
            # simulating fragment around read template
            fragPos = self._simFragPos(readTempPos)            
            
            # parsing fragment from genome index
            fragSeq = genome.get_seq(*fragPos)            
            fragSeqLen = len(fragSeq)
            
            # checking to see if fragment acutally in range (accounting for edge cases)
            if tryCnt >= 1000:
                raise ValueError('Exceeded {} tries to find frag length in min-max range'.format(str(tryCnt)))                
            elif (fragSeqLen >= self.get_minFragSize()
                and fragSeqLen <= self.get_maxFragSize()):
                break
            else:
                tryCnt += 1
                continue
                
        return [fragPos[0], fragPos[1], fragPos[2], fragSeq]

        
    def _simFragPos(self, readTempPos):
        """Simulate the position of the fragment containing the read template.
        Location selected using a uniform distribution.
        Fragment length selected using user-provide distribution (set during initilization)

        Args:
        readTempPos -- list: [scaffoldID,readStart,readTemplateEnd]
        Return:
        [scafID,fragStart,fragEnd]
        """
        readTempSize = readTempPos[2] - readTempPos[1] + 1
        
        # nfrag size
        tryCnt = 0
        while 1:
            fragSize = int(self.fld(size=1))
            #print 'fragSize: {}'.format(str(fragSize))
            if tryCnt >= 1000:
                raise ValueError('Exceeded {} tries to find frag length in min-max range'.format(str(tryCnt)))
            elif fragSize >= self.get_minFragSize() and fragSize <= self.get_maxFragSize():
                break
            else:
                tryCnt += 1
                continue
                

        # frag start (position upstream from read template)        
        randPos = int(pymcDist.runiform(0, fragSize - readTempSize, size=1))
        fragStart = readTempPos[2] - randPos
        fragEnd = fragStart + fragSize - 1
        
        return readTempPos[0], fragStart, fragEnd
                

    def _set_user_dist(self, userDist):
        """Setting user defined distribution. Using pymc distribution functions.
        Args:
        userDist -- User defined distribution with moment info. Example: ['normal',10,1]
        Return:
        function -- partial pymc function with moment values provided
        """
        userDist[0] = str(userDist[0]).lower()
        userDist[1:] = [int(x) for x in userDist[1:]]
        
        if userDist[0] == 'normal':
            assert len(userDist) >= 3, 'mu and tau must be provided for "normal" distribution'
            return functools.partial(np.random.normal, loc=userDist[1], scale=userDist[2])
        elif userDist[0] == 'uniform':
            assert len(userDist) >= 3, 'lower,upper must be provided for "uniform" distribution' 
            return functools.partial(np.random.uniform, low=userDist[1], high=userDist[2])
        elif userDist[0] == 'skewed-normal':
            assert len(userDist) >= 4, 'mu,tau,alpha must be provided for "skewed-normal" distribution'
            return functools.partial(skew_normal.rnd_skewnormal, location=userDist[1], scale=userDist[2], shape=userDist[3])
        else:
            raise ValueError('Distribution "{}" is not supported.'.format(userDist[0]))
                        
        
    def _get_randPos_amp(self, genome):
        """Getting the genomic position of a randomly selected amplicon.
        Args:
        genome -- genome-like object
        Return:
        [start,end] -- read template start-end
        """
        nAmps = genome.get_nAmplicons()
        if nAmps == 1:
            i = 0
        else:
            i = pymcDist.rdiscrete_uniform(0, nAmps - 1)
        row = genome.get_MFEprimerRes().iloc[i]
        row = row[['HitID','BindingStart','BindingStop']]
        return [x for x in row]
        

    def _get_randPos_shotgun(self, genome):
        """Randomly selecting start-end position of read template.
        start selected from a uniform distribution.
        end seleded based on read template length distribution.
        Args:
        genome -- genome-like object
        Return:
        [scafName,start,end] -- scaffold and start-end of read template
        """
        assert hasattr(self, 'rtl'), '"rtl" attribute required'
        
        # randomly selecting a genome scaffold
        scafName = random.choice(genome.fastaIdx.keys())
        scaf = genome.fastaIdx[scafName]

        # start
        start = pymcDist.rdiscrete_uniform(scaf.start, scaf.stop)
        
        # end        
        mpSize = self.rtl(size=1)
        assert mpSize > 0, 'Read template size must be > 0'

        # scaffoldName, readTempStart, readTempEnd
        return scafName, start, start + int(mpSize) -1

        
    def get_minFragSize(self):
        return self.flr[0]

    def get_maxFragSize(self):
        return self.flr[1]
                                    
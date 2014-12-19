"""Class for simulating fragments from genomes"""

# import
## batteries
import sys, os
import subprocess
from StringIO import StringIO
import itertools

## 3rd party
import pyfasta
#import pydna
from Bio import SeqIO
import pandas as pd
from intervaltree import Interval, IntervalTree
import numpy as np

## appliation
import Utils
from SimFrags import SimFrags


            
class Genome(object):
    """subclass of genomes: 1 genome entry"""

    def __init__(self, inFile, taxonName, primerFile=None):

        """
        Args:
        inFile -- name of the genome sequence file
        taxonName -- name of taxon
        primerFile -- file name of primers in fasta format
        """

        self.fileName = inFile
        self.taxonName = taxonName
        self.primerFile = primerFile

        # checking that the genome files exist
        Utils.checkExists(self.fileName)
        
        # create fasta index
        self.fastaIdx = pyfasta.Fasta(self.get_fileName())
            
        
    def callMFEprimer(self, rtr, MFEprimerExe='./MFEprimer.py'):
        """Calling MFEprimer to get amplicons.
        Args:
        rtr -- read template size range (min,max)
        MFEprimerExe -- string with path of MFEprimer.py
        """
        
        # Input check
        try:
            (minTemplateLen,maxTemplateLen) = rtr[:3]
        except IndexError:
            raise IndexError('"rtr" must contain (min,max)')

            
        # system call
        cmd = '{} -i {} -d {} --tab --size_start {} --size_stop {} --ppc 70 --tm_start 50 --tm_stop 70'
        cmd = cmd.format(MFEprimerExe,
                         self.get_primerFile(),
                         self.get_fileName(),
                         minTemplateLen,
                         maxTemplateLen)
        ret = subprocess.check_output(cmd, shell=True)
        
        # load results as dataframe
        self.MFEprimerRes = pd.read_csv(StringIO(ret), sep='\t')

        # editing values
        self.MFEprimerRes.BindingStart.astype(int)
        self.MFEprimerRes.BindingStop.astype(int)        

        
    def filterOverlaps(self, overlapPercCutoff=70):
        """Filtering out amplicons that substantially overlap.
        The amplicon with the highest PPC with be kept.
        The MFEprimerRes attribute must be set.
        
        Args:
        overlapPercCutoff -- percent of overlap to consider 'substantially' overlapping
        Attribute edits:
        MFEprimer table filtered of overlaps
        """
        if self.get_MFEprimerRes() is None:  
            raise AttributeError('genome object does not have MFEprimerRes attribute. Run MFEprimer() first')

            
        # making interval tree
        tree = IntervalTree()
        
        # loading intervals
        for count,row in self.get_MFEprimerRes().iterrows():
            # sanity check for + strand
            if row['BindingStart'] > row['BindingStop']:
                raise TypeError('MFEprimer binding start-stop is not + strand')
            tree.addi(row['BindingStart'], row['BindingStop'],
                      [count, row['PPC'], row['Size']])
                        
        # finding all that substantially overlap; keeping one with > PPC
        tree2 = tree.copy()
        for iv1 in tree.iter():
            # skipping if already removed from tree2
            if not iv1 in tree2:
                continue
                
            overlaps = tree.search(iv1.begin, iv1.end)

            # skipping those that poorly overlap 
            lowOverlap = set()
            for iv2 in overlaps:
                if iv1.overlaps(iv2):
                    percOverlaps = self._calcPercOverlap(iv1, iv2)
                    if percOverlaps[0] < overlapPercCutoff:
                        lowOverlap.add(iv2)
            overlaps = overlaps - lowOverlap  # just list of substantially overlapping

            # skipping those that have been already removed
            prevRm = set([x for x in overlaps if x not in tree2])
            overlaps = overlaps - prevRm
            
            # removing all substantially overlapping intervals with lower PPC
            if len(overlaps) > 1:
                overlaps = sorted(overlaps, key=lambda x: x.data[1], reverse=True)
                for o in overlaps[1:]:
                    if o in tree2:
                        tree2.remove(o)
            else:
                pass
 
        # selecting columns
        iv_idx = [x.data[0] for x in tree2.iter()]
        self.MFEprimerRes = self.MFEprimerRes.iloc[iv_idx]                
                      

    @staticmethod
    def _calcPercOverlap(iv1, iv2):
        """Calculating overlap between intervals (iv1, iv2)
        Return:
        tuple of percent overlap relative to iv1, relative to iv2
        """
        if not iv1.overlaps(iv2):
            return [0.0, 0.0]
        elif iv1 == iv2:
            return [100.0, 100.0]
        else:
            overlapSize = len(set(range(iv1.begin,iv1.end +1)) &
                              set(range(iv2.begin,iv2.end +2)))
            overlapSize = float(overlapSize)
            iv1PercOverlap = (iv1.end - iv1.begin + 1) / overlapSize * 100.0
            iv2PercOverlap = (iv2.end - iv2.begin + 1) / overlapSize * 100.0
            return [iv1PercOverlap, iv2PercOverlap]


    @staticmethod
    def calcGC(seq):
        """Calculating GC content of a sequence string. Return as percent G+C.
        Only unambiguous nucleotide codes will be counted.
        """
        seq = str(seq).upper()
        aCount = seq.count('A')
        tCount = seq.count('T')        
        cCount = seq.count('C')
        gCount = seq.count('G')
        
        return float(cCount + gCount) / float(aCount + tCount + cCount + gCount) * 100.0
            
                    
    # getters/setters/iters
    def get_fileName(self, rmPath=False):
        if rmPath is True:
            return os.path.split(self.fileName)[1]
        else:
            return self.fileName
            
    def get_taxonName(self):
        return self.taxonName
        
    def get_nFrags(self):
        return self.nFrags
        
    def get_minTemplateRange(self):
        return self.rtr[0]
        
    def get_maxTemplateRange(self):
        return self.rtr[1]
        
    def get_nAmplicons(self):
        try:
            return self.MFEprimerRes.shape[0]
        except ValueError:
            return None

    def get_primerFile(self):
        try:
            return self.primerFile
        except AttributeError:
            return None
                    
    def get_MFEprimerRes(self):
        """Getting the results of MFEprimer. Returns a dataframe"""
        try:
            return self.MFEprimerRes
        except AttributeError:
            return None
            
    def iter_seqRecs(self):
        """Iterate over sequence records"""
        with open(self.fileName) as inF:
            for rec in SeqIO.parse(inF, self.get_fileFormat()):
                yield rec
                
    def get_seq(self, scaffold, start, end):
        """Getting sequence from genome. 0-indexing.
        Args:
        scaffold -- scaffold id
        start -- sequence start position
        end -- sequence end position
        """
        try:
            return str(self.fastaIdx[scaffold][start:end])
        except AttributeError:
            raise AttributeError('No fastaIdx attribute for genome object')
        except KeyError:
            raise KeyError('ScaffoldID "{}" not found in genome fasta index'.format(scaffold))
    
"""Utility scripts for application"""

# import
## batteries
import os,sys
from pprint import pprint


def checkExists(f):
    """
    Args:
    f -- file name
    """
    if not os.path.isfile(f):
        raise IOError('"{}" not found. Did you provide the full PATH?'.format(f))


def parseGenomeList(inFile, filePath=None):
    """Parsing the genome list file
    Args:
    inFile -- genome list file name
    filePath -- abs path to genome sequence files
    Return:
    dict -- {genomeFile : taxonName}
    """
    # parse file as list
    genomeList = []
    with open(inFile, 'rb') as inF:
        for line in inF:
            row = line.rstrip().split('\t')
            
            if len(row) < 2:
                raise ValueError('Need format: "taxonName<tab>fileName";'
                                 'for row: "{}"'.format(row))
            else:
                (taxonName,fileName) = row[:2]
                
            # path to genome file
            if filePath is not None:
                fileName = os.path.join(filePath, fileName)

            # checking for file existence
            checkExists(fileName)
                
            #genomeList[fileName] = taxonName
            genomeList.append((fileName,taxonName))
                
    return genomeList


def describe_builtin(obj):
    """ Describe a builtin function """

    #wi('+Built-in Function: %s' % obj.__name__)
    # Built-in functions cannot be inspected by
    # inspect.getargspec. We have to try and parse
    # the __doc__ attribute of the function.
    docstr = obj.__doc__
    args = ''
    
    if docstr:
        items = docstr.split('\n')
#        pprint(items)
        if items:
            func_descr = items[0]
            s = func_descr.replace(obj.__name__,'')
            idx1 = s.find('(')
            idx2 = s.find(')',idx1)
            if idx1 != -1 and idx2 != -1 and (idx2>idx1+1):
                args = s[idx1+1:idx2]
                #wi('\t-Method Arguments:', args)
                for arg in args:
                    yield arg
                
    if args=='':
        yield None

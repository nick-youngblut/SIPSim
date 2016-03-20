"""Utility scripts for application"""

# import
## batteries
import os,sys
import re
import time
import platform
import subprocess
from pprint import pprint
from itertools import chain
from functools import partial
import dill
import random

## 3rd party
import numpy as np
import pandas as pd
import scipy.stats as ss


def get_os():
    """Get operating system; only works for unix-like machines"""
    OS = platform.uname()[0]
    if OS == 'Linux':
        OS = 'linux'
    elif OS == 'Darwin':
        OS = 'mac'
    else:
        sys.stderr.write('OS: "{}" not supported\n'.format(OS))

    return OS


def is_file(fileName):
    """Does file exist?"""
    if os.path.isfile(fileName) is False:
        raise IOError('"{}" does not exist'.format(fileName))

        
def sys_call(cmd, quiet=False):
    """System call of command.
    
    Parameters
    ----------
    cmd : str
        The command to run as a system call
    quiet : bool
        Suppress the system command output

    Returns
    -------
    output : str
       system call output
    err : str
       system call error
    """
    try:
        if quiet:
            DEVNULL = open(os.devnull, 'w')
            proc = subprocess.Popen([cmd], shell=True, stdout=DEVNULL)
        else:
            proc = subprocess.Popen([cmd], shell=True)
    except subprocess.CalledProcessError:
        pass # handle errors in the called executable
    except OSError:
        raise OSError('No executable for command: "{}"\n'.format(cmd))

    output, err = proc.communicate()


def load_kde(fileName):
    """Load a pickled dict {taxon:kde_object} file.
    
    Parameters
    ----------
    fileName : str
        name of pickled file ('-' if from STDIN)

    Returns
    -------
    dict : {taxon_name:kde_object}
    """
    try:
        if fileName == '-':
            kde = dill.load(sys.stdin)
        else:
            with open(fileName, 'rb') as inFH:
                kde = dill.load(inFH)
    except dill.UnpicklingError:
        msg = 'Cannot unpickle "{}"'
        raise dill.UnpicklingError, msg.format(fileName)

    return kde


def KDE_type(KDE_obj):
    """Determining the KDE object structure. Possible type:
    1) [taxon,kde]
    2) {taxon:kde}
    3) {libID:{taxon:kde}}
    
    Args
    ----
    KDE_obj : interable
        Some sort of KDE object used by SIPSim
    
    Returns
    -------
    int : number corresponding to KDE object type
    """
    kde_type = None
    try:            
        for x,y in KDE_obj.items():
            try: 
                for xx,yy in y.items():  # {libID:{taxon:kde}}
                    kde_type = 3
                    break
            except AttributeError:       # {taxon:kde}                
                kde_type = 2
                break
    except AttributeError:  # [taxon,kde]
        for x in KDE_obj:   
            try:
                for y in x:
                    kde_type = 1
                    break
            except TypeError:
                raise TypeError, 'KDE object type not recognized'
            break
    return kde_type
    

def checkExists(f):
    """ Check that the file `f` exists."""
    if not os.path.isfile(f):
        msg = '"{}" not found. Did you provide the full PATH?'
        raise IOError(msg.format(f))


def checkEmpty(f):
    """ Check that the file `f` is not empty"""
    if os.stat(f).st_size == 0:
        msg = '"{}" is empty!'
        raise IOError(msg.format(f))


def parseGenomeList(inFile, filePath=None, check_exists=True):
    """Parsing the genome list file.

    Parameters
    ----------
    inFile : str
        file name of genome list file
    filePath : str
        The absolute path to genome sequence files. a
    check_exists : bool
        Check if genome sequence files exist
    
    Returns:
    list : [[taxonName, genomeFile], ...]
    """
    # parse file as list
    genomeList = []
    with open(inFile, 'rb') as inF:
        for line in inF:
            row = line.rstrip().split('\t')
            
            if row[0] == '' or row[1] == '':
                raise IOError, "Necessary row value is empty!"

            if len(row) < 2:
                raise ValueError('Need format: "taxonName<tab>fileName";'
                                 'for row: "{}"'.format(row))
            else:
                (taxonName,fileName) = row[:2]
                
            # path to genome file
            if filePath is not None:
                fileName = os.path.join(filePath, fileName)

            # checking for file existence
            if check_exists:
                checkExists(fileName)
                
            #genomeList[fileName] = taxonName
            genomeList.append((taxonName,fileName))
                
    return genomeList


def describe_builtin(obj):
    """ Describe a builtin function if obj.__doc__
    available.

    Parameters
    ----------
    obj : python object
    
    Returns
    -------
    iterator : builtin args
    """
    #wi('+Built-in Function: %s' % obj.__name__)
    # Built-in functions cannot be inspected by
    # inspect.getargspec. We have to try and parse
    # the __doc__ attribute of the function.
    docstr = obj.__doc__
    args = ''
    
    if docstr:
        items = docstr.split('\n')
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


def parseKeyValueString(x):
    """Parse a string in format: 'key1:value1,key2:value2,keyN:valueN'.
    Values assumed to be numeric.
    
    Parameters
    ----------
    x : string
        Required format: 'key:value,key:value,...' or 'key=value,key=value,...'

    Returns
    -------
    dict : {key:value, ...}
    """
    if x is None or x == 'None':
        return {}
    x = x.replace(' ','')
    l = re.split('[=:,]', x)
    return {k.lower():float(v) for k,v in zip(l[::2],l[1::2])}


def random_walk_var_step(x, max_walk):
    """Shuffle the order of a list based on a random walk along list value
    ranks, where walk step size is randomly selected to be between 0 and
    max_walk for each step. 
    This produces a ordering with a certain level of autocorrelation
    that is set by the max walk step size. The larger the max walk step size,
    the less autocorrelation that will be prevalent. 
    
    Parameters
    ----------
    x : list
        Rank values
    max_walk : int
        The max distance in rank for the walk step.

    Returns
    -------
    list : a reordered list of values
    """
    # x as list
    try:
        x = list(x)
    except TypeError:
        msg = 'x must be a list-like object'
        raise TypeError, msg
        
    x_len = len(x)
    
    # ranks of values
    ## (index, value, value_rank)
    ranks = zip(range(x_len), x, ss.rankdata(x))
 
    # starting range
    start_i = random.randrange(0,x_len,1)
    cur_rank = ranks[start_i]

    # filtering cur_rank from ranks 
    ranks = [x for x in ranks if x[0] != cur_rank[0]]
    
    # moving through ranks
    x_new_order = [cur_rank[1]]
    for i in xrange(x_len-1):       
        # select max walk distance
        if max_walk > 1:
            max_range = random.randrange(1,max_walk)
        else:
            max_range = 1
    
        # filter to just ranks w/in rank distance
        filt_ranks = [x for x in ranks 
                      if abs(x[2] - cur_rank[2]) <= max_range]
    
        # selecting randomly from filtered ranks
        cur_rank = random.sample(filt_ranks, k=1)[0]
        
        # adding value to new list
        x_new_order.append(cur_rank[1])
        
        # filtering cur_rank from ranks 
        ranks = [x for x in ranks if x[0] != cur_rank[0]]

        # re-ranking remaining values
        rank_vals = [x[2] for x in ranks]
        new_ranks = ss.rankdata(rank_vals)
        for i,v in enumerate(new_ranks):
            ranks[i] = (ranks[i][0], ranks[i][1], v)
                    
    return x_new_order


def part_dist_func(dist, dist_params):
    """Creating a numpy.random distribution function with the
    distribution parameters already set.
 
    Parameters
    ----------
    dist : str
        name of numpy.random distribution function
    dist_params : dict
        numpy.random distribution function parameters.
        Example: {low:1, high:10}

    Returns
    -------
    function : numpy.random function with set parameters
    """
    ## get numpy function
    try:
        dist_func =  getattr(np.random, dist)
    except AttributeError:
        raise AttributeError('Distribution "{}" not supported\n'.format(dist))

    # if function should return one constant value
    try:
        if dist_params['low'] == dist_params['high']:
            return lambda size: [dist_params['low']] * size
    except KeyError:
        pass

    # else making partial function
    try:
        part_dist_func = partial(dist_func, **dist_params)
        part_dist_func(size=1)
    except TypeError:
        params = ','.join([str(x) + ':' + str(y)  for x,y 
                           in dist_params.items()])
        msg = 'Params "{}" do not work with distribution "{}"\n'
        raise TypeError(msg.format(params, dist))  

    return part_dist_func

        

class Status(object):
    """Simple custom logging information"""
    def __init__(self, quiet=False):
        self.quiet = quiet
        self.msgs = {'kde':'GC/fragment_length KDE sampled',
                     'diffusion':'diffusion added to BD values',
                     'incorp':'isotope incorporation added to BD values',
                     'bin':'binned BD values',
                     'final':'taxon finished',
                     'zero':'NOTE: taxon has an abundance of 0'}                

    def msg(self, msgKey, startTime=None):
        """Writing formatted status message to STDERR.
        
        Parameters
        ----------
        msgKey : str
            dict key for writing which status message
        startTime : bool
            used to measure how much time has elapsed        
        """    
        nowTime = time.time()
        if startTime is not None:
            timeDiff = '{0:.1f}'.format(nowTime - startTime)
        else:
            timeDiff = '0.0'

        if not self.quiet:
            try:
                x = '     Elapsed: {0:>7} sec => {1}\n'
                sys.stderr.write(x.format(timeDiff, self.msgs[msgKey.lower()]))
            except KeyError:
                s = 'Cannot find status for message key "{}"'
                raise KeyError(s.format(msgKey))

        return nowTime

    def status_quiet(self):
        self.quiet = True

    def status_loud(self):
        self.quiet = False

    def get_msgKeys(self):
        """Get all possible message keys.
        """
        return self.msgs.keys()
        
        
class _table(object):
    """Template class for reading in SIPSim tables.
    Tables are just pandas.DataFrame objects.
    """    
    def __init__(self, df, filename):
        """
        Parameters
        ----------
        df : pandas dataframe object
        filename : str
             name of table file
        """
        self.df = df
        self.tableFileName = filename
        
        # library as string
        try:
            self.df['library'] = self.df['library'].astype(str)
        except KeyError:
            self.wide2long()
            try:
                self.df['library'] = self.df['library'].astype(str)
            except KeyError:
                msg = '"library" column not found in table: "{}"!'
                raise KeyError(msg.format(filename))
            

    # reshaping table
    def wide2long(self, sep='__'):
        """Convert table from wide to long format.
        
        Parameters
        ----------
        sep : str
            used to split column names
        """
        self.df = pd.melt(self.df, id_vars=['taxon'])
        new_cols = self.df['variable'].str.split(sep).apply(pd.Series, 1)
        d = {'value':'count',0:'library',1:'fraction'}
        self.df = self.df.join(new_cols)\
                         .rename(columns=d)\
                         .drop('variable',1)

        try:
            l = ['library','fractions','taxon','count']
            self.df = self.df\
                          .reindex_axis(l, axis=1)\
                          .sort_values(by=['taxon', 'fraction', 'library'])
        except KeyError:
            pass


    def long2wide(self, values, index, columns):
        """Convert table from long to wide format.
        
        Parameters
        ----------
        values : list
            values in pivot table
        index : list?
            index in pivot table
        columns : list
            columns in pivot table
        """
        self.df = pd.pivot_table(self.df, values=values, index=index,
                                 columns=columns, fill_value=0)
        self.df.columns =  ['__'.join(x) for x in self.df.columns.tolist()]

            
    # write table
    def to_csv(self, *args, **kwargs):
        self.df.to_csv(*args, **kwargs)

    def write_table(self, *args, **kwargs):
        self.to_csv(*args, **kwargs)

            
    # load from csv
    @classmethod
    def from_csv(cls, filename, **kwargs):
        """Read in table file to a pandas dataframe.

        Parameters
        ----------
        filename : str
            Table file name
        kwargs : dict
            passed to pandas.read_csv

        Returns
        -------
        pandas.DataFrame subclass 
        """
        df = pd.read_csv(filename, **kwargs)
        return cls(df, filename)

        
    # get/set/iter
    def iter_uniqueColumnValues(self, columnID):
        """General iteration of unique column values.

        Parameters
        ----------
        str : ID of column in table
        """
        try:
            for l in self.df[columnID].unique():
                yield l
        except KeyError:
            raise KeyError('Column "{}" not found'.format(columnID))
            
    def iter_libraries(self):
        """iterate through all unique library IDs."""
        for libID in self.iter_uniqueColumnValues('library'):
            yield libID
                
    def iter_taxa(self, libID=None):
        """Iterate through all unique taxon names."""
        col_name = None        
        try:
            self.df['taxon_name']
            col_name = 'taxon_name'
        except KeyError:
            try:
                self.df['taxon']
                col_name = 'taxon'
            except KeyError:
                raise KeyError('Neither "taxon_name" nor "taxon" is a column')
        if libID is None:            
            for taxon_name in self.iter_uniqueColumnValues(col_name):
                yield taxon_name
        else:
            df_lib = self.df.loc[self.df['library'] == libID]
            for taxon_name in df_lib[col_name].unique():
                yield taxon_name
                
    def iter_taxonRowsInLib(self, libID):
        """Iterate through all subset dataframes containing just 1 taxon.
        
        Parameters
        ----------
        libID : str
            library ID        
        """
        df_lib = self.df.loc[self.df['library'] == libID]
        for taxon_name in df_lib['taxon_name'].unique():
            yield (taxon_name, df_lib.loc[df_lib['taxon_name'] == taxon_name])
                
    def __repr__(self):
        return self.df.__repr__()
        

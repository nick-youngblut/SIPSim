#!/usr/bin/env python

#--- Option parsing ---#
"""
fragment_parse: parsing out fragment objects for certain genomes

Usage:
  fragment_parse [options] <fragments> <genome_list>
  fragment_parse [options] --random=<rn> <fragment_kde>
  fragment_parse -h | --help
  fragment_parse --version

Options:
  <fragments>     Simulate fragment file.
  <genome_list>   Table listing genomes to parse.
  --name=<n>      Column in genomes table that contains names of genome
                  fragments to parse out.
                  [default: genomeID]
  --rename=<r>    Column in genomes table to rename parsed genome fragments.
                  [default: None]
  --cluster=<c>   Column that designates matching clusters for selecting 
                  target genomes. 
                  [default: None]
  --NA-random     For NA values in --name column, parse out random 
                  fragments and rename by --rename column (must provide 
                  --rename).
  --random=<rn>   Parse out <rn> randomly selected genome fragments.
                  Sampling with replacement of genome fragments.
                  NOTE: genome_list file is ignored
                  [default: 0]
  --invert        Invert selection.
  --log=<l>       Name of log file.
                  [default: fragment_parse.log]
  --debug         Debug mode.
  --version       Show version.
  -h --help       Show this screen.

Description:
  Parsing out fragments from select genomes.
  Each '.' in the status output denotes 10 fragment objects parsed.
  
  invert
  ------
  Get the opposite selection of genome fragments from what is provided in the 
  genome list.

  random
  ------
  Get a random list of genome fragments. Fragments are drawn from list of 
  fragments without replacement, so you can create a list of genome fragmentss as 
  long as you want. The fragemnts will be named: "OTU.rand#"

  cluster
  -------
  All rows in genome_list that have 'NA' or equivalent in --cluster column 
  will be skipped.

  Output
  ------
  Fragments object written to STDOUT.
"""

# import
## batteries
from docopt import docopt
import sys,os
import copy
import dill
import numpy as np
import pandas as pd
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)


# functions
def parse_comm_by_cluster(comm, clust_col, inv=False):
    """Parsing community file.    
    Parameters
    ----------
    comm : pandas.dataframe
        community table
    clust_col : str
        column for selecting target OTUs
    inv : bool
        invert selection
    
    Returns
    -------
    comm : pandas.dataframe
    """    
    if(inv):
        comm = comm.loc[np.isnan(comm[clust_col])]        
    else:
        comm = comm.loc[comm[clust_col] > 0]
        comm = comm.dropna()
    # status
    msg = 'Genome list length (removed all NAs): {}\n'
    sys.stderr.write(msg.format(comm.shape[0]))

    return(comm)


def parse_target_frag_kdes(frag_kdes, comm, NA_random, log):
    """Parsing out genome fragment KDEs for target taxa.
    Parameters
    ----------
    frag_kdes : list 
        list of frag_kde objects
    comm : pandas.dataframe
        community table
    NA_random : bool
        NAs (non-targets) get a random fragment KDE
    log : buffer
        Log file buffer

    Returns
    -------
    frag_kdes
    """    
    # making an index for frag_kde list (genomeID : list_index)
    target_idx = {g[0]:i for i,g in enumerate(frag_kdes)}
    msg = 'Number of genomes with simulated fragments: {}\n'
    len_idx = len(target_idx)
    sys.stderr.write(msg.format(len_idx))
        
    # parsing out frag_kdes
    frag_kde_target = []
    stat_cnt = 0
    for count,row in comm.iterrows():
        stat_cnt += 1
        if stat_cnt % 10 == 0:
            sys.stderr.write('.')
        tg = row[0]   # genomeID by default
        # getting index location of KDE in KDE list 
        try:
            idx = target_idx[tg]
        except KeyError:
            if NA_random:
                # random selection of fragment_kde 
                idx = np.random.choice(len_idx,1)[0]
            else:
                idx = None
                msg = 'Cannot find "{}" in target list (check --name column)'
                print msg.format(tg)
        # copy KDE 
        try:
            frag_kde_target.append(copy.deepcopy(frag_kdes[idx]))
        except TypeError:
            continue
        # renaming KDE from rename column in comm file
        ## skipping if no rename ID exists
        try:
            idx2 = len(frag_kde_target)-1
            oldname = frag_kde_target[idx2][0]
            newname = row[1]   # --rename column in comm file
            frag_kde_target[idx2][0] = newname
            msg = 'Renamed fragment KDE: old={}, new={}\n'
            log.write(msg.format(oldname, newname))
        except IndexError:
            continue
            
    sys.stderr.write('\n')
                
    # check that all KDEs were found
    x_len = len(frag_kde_target)
    y_len = comm.shape[0]
    if(x_len != y_len):
        msg = 'WARNING: No matching genomes for {} OTUs!\n'
        sys.stderr.write(msg.format(y_len - x_len))

    return(frag_kde_target)   


def parse_nontarget_frag_kdes(frag_kdes, comm, log): 
    """Parsing out genome fragment KDEs for non-target taxa.
    Parameters
    ----------
    frag_kdes : list 
        list of frag_kde objects
    comm : pandas.dataframe
        community table
    log : buffer
        log file buffer
    
    Returns
    -------
    frag_kdes
    """
    richness_needed = comm.shape[0]
    if richness_needed <= 0:
        msg = 'Nothing in the parse list!\n'
        sys.stderr.write(msg)
        sys.exit()

    # making an index of non-target genome fragments
    ## will select randomly from index to select genome fragments
    nt = set(comm.iloc[:,0])
    frag_nonTarget_idx = [i for i,x in enumerate(frag_kdes) \
                             if x[0] not in nt]
    pool_size = len(frag_nonTarget_idx)
    msg = 'Non-target fragment pool size: {}\n'
    sys.stderr.write(msg.format(pool_size))
        
    # parsing amp frag
    ## random choice index
    r_idx = range(len(frag_nonTarget_idx))
    r_idx = np.random.choice(r_idx, richness_needed)    
    ## parsing out randomly selected fragments for each non-target OTU
    stat_cnt = 0
    frag_kde_rand = []
    for i in r_idx:
        stat_cnt += 1
        if stat_cnt % 10 == 0:
            sys.stderr.write('.')
        ii = frag_nonTarget_idx[i]
        frag_kde_rand.append(copy.deepcopy(frag_kdes[ii]))
    sys.stderr.write('\n')
            
    # renaming randomly selected KDEs by random OTU_IDs
    try:
        for i in range(richness_needed):
            oldname = frag_kde_rand[i][0]
            newname = comm.iloc[i,1]
            frag_kde_rand[i][0] = newname
            msg = 'Renamed fragment KDE: old={}, new={}\n'
            log.write(msg.format(oldname, newname))
    except IndexError:
        pass
        
    # checking
    x_len = len(frag_kde_rand)
    if(x_len != richness_needed):
        msg = 'Richness != needed richess; non-target richness = {}'
        sys.stderr.write(msg.format(x_len))       
            
    return frag_kde_rand


def parse_random_frag_kdes(frag_kdes, n, log):         
    """Sampling genome fragment KDEs from random genomes. 
    Subsampling at genome level; subsampling w/out replacement.
    Parameters
    ----------
    frag_kdes : list 
        list of frag_kde objects
    n : int
        Number of subsamples
    log : buffer
        Log file 

    Returns
    -------
    frag_kdes
    """
    # parsing amp frag
    ## random choice index
    r_idx = range(len(frag_kdes))
    r_idx = np.random.choice(r_idx, n)
    ## parsing out randomly selected fragments
    stat_cnt = 0
    frag_kde_rand = []
    for i in r_idx:
        stat_cnt += 1
        if stat_cnt % 10 == 0:
            sys.stderr.write('.')
        frag_kde_rand.append(copy.deepcopy(frag_kdes[i]))
    sys.stderr.write('\n')
            
    # renaming randomly selected KDEs by random OTU_IDs
    for i in range(richness_needed):
        oldname = frag_kde_rand[i][0]
        newname = 'OTU.rand{}'.format(i) 
        frag_kde_rand[i][0] = newname
        msg = 'Renamed fragment KDE: old={}, new={}\n'
        log.write(msg.format(oldname, newname))

    return frag_kde_rand
            

# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')
    args['--random'] = int(args['--random'])    

    # open log file
    logFH = open(args['--log'], 'wb') 

    # loading comm file
    if(args['--random'] <= 0):
        df_comm = pd.read_csv(args['<genome_list>'], sep='\t', na_values=['NA'])
        ## parsing rows of comm file
        if(args['--cluster'] != 'None'):
            df_comm = parse_comm_by_cluster(df_comm, 
                                            args['--cluster'], 
                                            inv=args['--invert'])
        ## parsing columns of comm file
        cols2parse = [args['--name']]
        if(args['--rename'] != 'None'):
            cols2parse += [args['--rename']]
        df_comm = df_comm.loc[:,cols2parse]


    # loading fragments
    sys.stderr.write('Loading fragments object...\n')
    frag_kdes = []
    with open(args['<fragments>'], 'rb') as iFH:
        frag_kdes = dill.load(iFH)

    # parsing fragments
    sys.stderr.write('Parsing fragments...\n')
    if(args['--random'] > 0):
        frag_kdes = parse_rand_frag_kdes(frag_kdes, args['--random'], logFH)
    elif(args['--invert']):
        frag_kdes = parse_nontarget_frag_kdes(frag_kdes, df_comm, logFH)
    else:
        frag_kdes = parse_target_frag_kdes(frag_kdes, df_comm, 
                                           args['--NA-random'], logFH)
        
    # writing parsed fragment file
    dill.dump(frag_kdes, sys.stdout)

    # closing log file
    logFH.close()
    msg = 'Log file written: {}\n'
    sys.stderr.write(msg.format(args['--log']))

#!/usr/bin/env python

"""
fragment_KDE_cat: concatenating 2 fragment_kde objects

Usage:
  fragment_KDE_cat [options] <fragment_kde1> <fragment_kde2>
  fragment_KDE_cat -h | --help
  fragment_KDE_cat --version

Options:
  <fragment_kde1>  Output from the fragment_kde subcommand.
  <fragment_kde2>  Output from the fragment_kde subcommand.
  --debug          Debug mode (turn off parallel processing).
  --version        Show version.
  -h --help        Show this screen.

Description:
  Concatenating fragment_kde objects produced by the fragment_kde subcommand.
  One list of fragment KDEs is appended to the other.

  Output
  ------
  Fragment KDE object written to STDOUT.
"""

# import
## batteries
from docopt import docopt
import sys,os
import dill

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)

    # load kde objects
    with open(args['<fragment_kde1>'], 'rb') as iFH1:
        kde1 = dill.load(iFH1)
    with open(args['<fragment_kde2>'], 'rb') as iFH2:
        kde2 = dill.load(iFH2)

    # joining and writing
    try:
        dill.dump(kde1 + kde2, sys.stdout)
    except TypeError:
        dill.dump(merge_two_dicts(kde1, kde2), sys.stdout)

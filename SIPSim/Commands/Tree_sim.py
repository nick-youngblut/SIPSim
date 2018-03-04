#!/usr/bin/env python

"""
tree_sim: Simulate a phylogeny for a specified set of taxa

Usage:
  tree_sim [options] list <genome_list>
  tree_sim [options] comm <comm_file>
  tree_sim -h | --help
  tree_sim --version

Options:
  <genome_list>          A file listing: taxonName<tab>genomeSeqFileName
  <comm_file>            Name of file produced by gradientComms subcommand.
  --birth_rate=<br>      Birth rate for birth-death model. 
                         [default: 1.0]
  --death_rate=<dr>      Death rate for birth-death model.
                         [default: 0.5]
  --birth_rate_sd=<bs>   Birth rate standard deviation for birth-death model. 
                         [default: 0.0]
  --death_rate_sd=<ds>   Death rate standard deviation for birth-death model. 
                         [default: 0.0]
  --star                 Create a star tree instead of a birth-death model tree.
  --outfmt=<f>           Output file format ('newick' or 'nexus').
                         [default: newick] 
  --version              Show version.
  --debug                Debug mode

Description:
  Simulate a tree based on a list of taxa provided either by a <genome_list>
  file or a <comm_file>.

  A simple birth-death model is used for tree simulation.

  Output
  ------
  A phylogeny is written to STDOUT.
"""

# import
## batteries
from docopt import docopt
import os, sys
## 3rd party
import dendropy
try:
    from dendropy.treesim import star_tree, birth_death    
except ImportError:
    from dendropy.simulate.treesim import star_tree
    from dendropy.model.birthdeath import birth_death_tree as birth_death
## application libraries
from SIPSim import Utils
from SIPSim.CommTable import CommTable


def main(args=None):    
    for param in ['birth_rate','death_rate','birth_rate_sd','death_rate_sd']:
        param = '--' + param
        args[param] = float(args[param])

    # loading taxon list 
    if args['<genome_list>'] is not None:
        taxa = Utils.parseGenomeList(args['<genome_list>'], check_exists=False)
        taxa = [x[0] for x in taxa]
    elif args['<comm_file>'] is not None:
        comm = CommTable.from_csv(args['<comm_file>'], sep='\t')
        taxa = comm.get_unique_taxon_names()

    # init dendropy taxon namespace
    taxa = dendropy.TaxonNamespace(taxa, label='taxa')

    # simulating tree
    if args['--star']:
        tree = star_tree(taxon_set=taxa)
    else:
        tree = birth_death(args['--birth_rate'],
                           args['--death_rate'],
                           birth_rate_sd=args['--birth_rate_sd'],
                           death_rate_sd=args['--death_rate_sd'],
                           num_extant_tips=len(taxa))
        
    # writing tree
    outfmt = args['--outfmt'].lower()
    psbl_fmts = ['newick','nexus']
    assert outfmt in psbl_fmts, 'output file format not recognized.' +\
        ' Possible formats: {}'.format(', '.join(psbl_fmts))
    tree.write_to_stream(sys.stdout, outfmt)
    
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)


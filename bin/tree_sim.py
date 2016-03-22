#!/usr/bin/env python

#--- Option parsing ---#
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
from dendropy.treesim import star_tree, birth_death
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
import Utils
from CommTable import CommTable


# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    
    for param in ['birth_rate','death_rate','birth_rate_sd','death_rate_sd']:
        param = '--' + param
        Uargs[param] = float(Uargs[param])

    # loading taxon list 
    if Uargs['<genome_list>'] is not None:
        taxa = Utils.parseGenomeList(Uargs['<genome_list>'], check_exists=False)
        taxa = [x[0] for x in taxa]
    elif Uargs['<comm_file>'] is not None:
        comm = CommTable.from_csv(Uargs['<comm_file>'], sep='\t')
        taxa = comm.get_unique_taxon_names()

    taxa = dendropy.TaxonSet(taxa)

    # simulating tree
    if Uargs['--star']:
        tree = star_tree(taxon_set=taxa)
    else:
        tree = birth_death(Uargs['--birth_rate'],
                           Uargs['--death_rate'],
                           birth_rate_sd=Uargs['--birth_rate_sd'],
                           death_rate_sd=Uargs['--death_rate_sd'],
                           taxon_set=taxa)
        
    # writing tree
    outfmt = Uargs['--outfmt'].lower()
    psbl_fmts = ['newick','nexus']
    assert outfmt in psbl_fmts, 'output file format not recognized.' +\
        ' Possible formats: {}'.format(', '.join(psbl_fmts))
    tree.write_to_stream(sys.stdout, outfmt)
    

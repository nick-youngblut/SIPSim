#!/usr/bin/env python

# import
## batteries
import os
import sys
## 3rd party
from docopt import docopt
#import Utils
## application
from SIPSim.Commands import BD_Shift
from SIPSim.Commands import Communities
from SIPSim.Commands import DBL
from SIPSim.Commands import DeltaBD
from SIPSim.Commands import Diffusion
from SIPSim.Commands import Fragment_KDE_cat
from SIPSim.Commands import Fragment_KDE
from SIPSim.Commands import Fragment_parse
from SIPSim.Commands import Fragments
from SIPSim.Commands import Genome_index
from SIPSim.Commands import Genome_rename
from SIPSim.Commands import Gradient_fractions
from SIPSim.Commands import HRSIP
from SIPSim.Commands import IncorpConfigExample
from SIPSim.Commands import Isotope_incorp
from SIPSim.Commands import KDE_bandwidth
from SIPSim.Commands import KDE_info 
from SIPSim.Commands import KDE_parse
from SIPSim.Commands import KDE_plot
from SIPSim.Commands import KDE_sample
from SIPSim.Commands import KDE_selectTaxa
from SIPSim.Commands import OTU_add_error
from SIPSim.Commands import OTU_PCR
from SIPSim.Commands import OTU_sampleData
from SIPSim.Commands import OTU_subsample
from SIPSim.Commands import OTU_sum
from SIPSim.Commands import OTU_table
from SIPSim.Commands import OTU_wideLong
from SIPSim.Commands import qSIP_atomExcess
from SIPSim.Commands import qSIP
from SIPSim.Commands import Tree_sim

def main(args=None):
    """Main entry point for application
    """
    if args is None:
        args = sys.argv[1:]
    
    docs = """
SIPSim: simulate isopycnic gradient fractionation of microbial community DNA

Usage:
  SIPSim <command> [<args>...]
  SIPSim -l | --list
  SIPSim -h | --help
  SIPSim --version

Options:
  -l --list     List subcommands.
  -h --help     Show this screen.
  --version     Show version.

Commands:
  Use the `list` option.
Description:
  Simulate how taxa would be distributed in isopycnic gradients as assessed by
  high throughput sequencing.
    """
    # arg parse
    args = docopt(docs,
                  version='0.1',
                  options_first=True)

    # dict of all subcommands
    cmds = {'BD_shift' : BD_Shift,
            'communities' : Communities,
            'DBL' : DBL,
            'deltaBD' : DeltaBD,
            'diffusion' : Diffusion,
            'fragment_KDE_cat' : Fragment_KDE_cat,
            'fragment_KDE' : Fragment_KDE,
            'fragment_parse' : Fragment_parse,
            'fragments' : Fragments,
            'genome_index' : Genome_index,
            'genome_rename' : Genome_rename,
            'gradient_fractions' : Gradient_fractions,
            'HRSIP' : HRSIP,
            'incorp_config_example' : IncorpConfigExample,
            'isotope_incorp' : Isotope_incorp,
            'KDE_bandwidth' : KDE_bandwidth,
            'KDE_info' : KDE_info,
            'KDE_parse' : KDE_parse,
            'KDE_plot' : KDE_plot,
            'KDE_sample' : KDE_sample,
            'KDE_select_taxa' : KDE_selectTaxa,
            'OTU_add_error' : OTU_add_error,
            'OTU_PCR' : OTU_PCR,
            'OTU_sample_data' : OTU_sampleData,
            'OTU_subsample' : OTU_subsample,
            'OTU_sum' : OTU_sum,
            'OTU_table' : OTU_table,
            'OTU_wide_long' : OTU_wideLong,
            'qSIP_atom_excess' : qSIP_atomExcess,
            'qSIP' : qSIP,
            'tree_sim' : Tree_sim}
    
    # list subcommands
    if args['--list']:
        cmd_list = '\n'.join(sorted(cmds.keys(), key=str.lower))
        print('#-- Commands --#')
        print(cmd_list)
        exit()

    # running subcommand
    try:
        cmds[args['<command>']].opt_parse(args['<args>'])
    except KeyError:
        msg = 'ERROR: command "{}" does not exist'
        print(msg.format(args['<command>']))
        exit()
        
    
if __name__ == '__main__':
    main()

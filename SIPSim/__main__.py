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
            'Difusion' : Diffusion,
            'fragment_KDE_cat' : Fragment_KDE_cat,
            'fragment_KDE' : Fragment_KDE,
            'fragment_parse' : Fragment_parse,
            'fragments' : Fragments,
            'genome_index' : Genome_index}
    
    # list subcommands
    if args['--list']:
        cmd_list = '\n'.join(sorted(cmds.keys()))
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

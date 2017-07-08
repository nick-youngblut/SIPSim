#!/usr/bin/env python

# import
## batteries
import os
import sys
## 3rd party
from docopt import docopt
#import Utils
## application
from SIPSim import BD_Shift


def main(args=None):
    """Main entry point for application
    """
    if args is None:
        args = sys.argv[1:]
    
    docs = """
SIPSim: simulate gradient fractionation of microbial community DNA

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
    cmds = {'BD_shift' : BD_Shift}
    
    # list subcommands
    if args['--list']:
        #Utils.list_subcmds(subcmdDir, '.py')
        cmd_list = ',\n'.join(cmds.keys())
        print('#-- Commands --#')
        print(cmd_list)
        exit()

    # running subcommand
    try:
        cmds[args['<command>']].opt_parse(args['<args>'])
    except KeyError:
        msg = 'Command "{}" does not exist'
        print(msg.format(args['<command>']))
        exit()
        
#    if args['<command>'] == 'BD_shift':
#        print(args['<args>'])
#        BD_Shift.opt_parse(args['<args>'])
        
    # calling subcommand scripts
    #scriptFile = os.path.join(subcmdDir, args['<command>'] + '.py')
    #if os.path.isfile(scriptFile):
    #    exit(call(['python', scriptFile] + args['<args>']))    
    #else:
    #    msg = '"{}" is not a SIPSim command. See "SIPSim -h".'
    #    exit(msg.format(args['<command>']))
    
if __name__ == '__main__':
    main()

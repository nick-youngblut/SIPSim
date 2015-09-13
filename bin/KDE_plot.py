#!/usr/bin/env python


#--- Option parsing ---#
"""
KDE_plot: make plots of each KDE (1D or 2D)

Usage:
  KDE_plot [options] <kde>
  KDE_plot -h | --help
  KDE_plot --version

Options:
  <kde>         Pickled KDE object
                ('-' if input from STDIN) 
  -o=<o>        Output file name.
                [Default: KDE.png]
  -n=<n>        Number of taxon KDEs to plot (0 = all plotted).
                [Default: 0]
  --nCol=<nc>   Number of subplot columns.
                [Default: 1]
  --xMin=<xm>   Minimum x-axis value ('' = min value in dataset).
                [Default: ]
  --xMax=<xM>   Maximum x-axis value ('' = max value in dataset).
                [Default: ]
  --yMin=<ym>   Minimum y-axis value ('' = min value in dataset).
                [Default: ]
  --yMax=<yM>   Maximum y-axis value ('' = max value in dataset).
                [Default: ]
  --xStep=<xs>  X dimension granularity.
                [Default: 0.0005]
  --yStep=<yx>  Y dimension granularity.
                [Default: 100]
  --xX=<xx>     X dimension figure size multiplier (ncol * x)
                [Default: 4]
  --yX=<yx>     Y dimension figure size multiplier (ncol * x)
                [Default: 3.5]
  --logY=<ly>   Base for y-axis log scaling ('' = no log scaling).
                [Default: ]
  -h --help     Show this screen.
  --version     Show version.
  --debug       Debug mode

Description:
  Plot each KDE (1D or 2D KDEs) in the the provided multi-KDE object. 

  Output:
  Image files written to `-o`
"""

# import
## batteries
from docopt import docopt
import sys,os
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)
# application
import Utils
import FigGen


# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')
    KDEs = Utils.load_kde(args['<kde>'])

    try:
        FigGen.KDE_ndims(KDEs)
    except AttributeError:
        outFile = os.path.splitext(args['-o'])
        msg = 'Processing library: "{}"\n'
        for lib,x in KDEs.items():
            sys.stderr.write(msg.format(lib))
            outName = ''.join([outFile[0], '_', str(lib), outFile[1]])
            FigGen.make_kde_fig(x, outName, 
                                n_subplot=args['-n'], 
                                ncol=args['--nCol'], 
                                xMin=args['--xMin'],
                                xMax=args['--xMax'],
                                yMin=args['--yMin'],
                                yMax=args['--yMax'],
                                xStep=args['--xStep'], 
                                yStep=args['--yStep'],
                                xX=args['--xX'],
                                yX=args['--yX'],
                                logY=args['--logY'])
    else:
        FigGen.make_kde_fig(KDEs, args['-o'], 
                            n_subplot=args['-n'], 
                            ncol=args['--nCol'], 
                            xMin=args['--xMin'],
                            xMax=args['--xMax'],
                            yMin=args['--yMin'],
                            yMax=args['--yMax'],
                            xStep=args['--xStep'], 
                            yStep=args['--yStep'],
                            xX=args['--xX'],
                            yX=args['--yX'],
                            logY=args['--logY'])

        

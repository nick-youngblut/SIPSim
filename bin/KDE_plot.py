#!/usr/bin/env python


#--- Option parsing ---#
"""
KDE_plot: make plots of each KDE (1D or 2D)

Usage:
  fragments [options] <kde>
  fragments -h | --help
  fragments --version

Options:
  <kde>         Pickled KDE object
                ('-' if input from STDIN) 
  -o=<o>        Output file name.
                [Default: KDE.png]
  -n=<n>        Number of taxon KDEs to plot (0 = all plotted).
                [Default: 0]
  --nCol=<nc>   Number of subplot columns.
                [Default: 3]
  --xStep=<xs>  X dimension granularity.
                [Default: 0.0005]
  --yStep=<yx>  Y dimension granularity.
                [Default: 100]
  --xX=<xx>     X dimension figure size multiplier (ncol * x)
                [Default: 4]
  --yX=<yx>     Y dimension figure size multiplier (ncol * x)
                [Default: 3.5]
  -h --help     Show this screen.
  --version     Show version.
  --debug       Debug mode

Description:
  Plot each KDE (1D or 2D KDEs) in the the provided kde object. 
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


    FigGen.make_kde_fig(KDEs, args['-o'], 
                        n_subplot=args['-n'], 
                        ncol=args['--nCol'], 
                        xStep=args['--xStep'], 
                        yStep=args['--yStep'],
                        xX=args['--xX'],
                        yX=args['--yX'])


    # checking n-dims    
    # ndims = FigGen.KDE_ndims(KDEs)
    # if len(ndims) > 1:
    #     msg = 'KDEs do not all have the same number of dimensions'
    # else:
    #     ndims = list(ndims)[0]

    # # plotting based on ndims
    # if ndims == 1:
    #     # making histograms of the 1d KDEs
    #     FigGen.make_kde_fig(KDEs, args['-o'], 
    #                              n_subplot=args['-n'], 
    #                              ncol=args['--nCol'], 
    #                              xStep=args['--xStep'], 
    #                              yStep=args['--yStep'],
    #                              xX=args['--xX'],
    #                              yX=args['--yX'])
    # elif ndims == 2:
    #     # making heatmaps of the 2d KDEs
    #     FigGen.make_2d_kde_plots(KDEs, args['-o'], 
    #                              n_subplot=args['-n'], 
    #                              ncol=args['--nCol'], 
    #                              xStep=args['--xStep'], 
    #                              yStep=args['--yStep'],
    #                              xX=args['--xX'],
    #                              yX=args['--yX'])
    # else:
    #     msg = 'KDE ndims > 2; cannot create figure'
    #     raise ValueError, msg
    


   
        

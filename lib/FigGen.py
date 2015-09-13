import sys,os
import math
import itertools
from collections import Counter
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def KDE_ndims(KDEs):
    """Get the number of dimesions of each KDE in the dict of KDEs
    
    Parameters
    ----------
    KDEs : dict
        dict of KDE objects
    
    Returns
    -------
    ndims : set       
    """
    c = Counter([KDE.d for KDE in KDEs.values() if KDE is not None])
    return set(c.keys())


def format_subplot_axes(subplot, xStart, xEnd, yStart=None, yEnd=None):
    """Formatting subplot axes to look pretty.
    
    Parameters
    ----------
    subplot : matplotlib subplot object
    xStart : float
        x-axis start
    xEnd : float
        x-axis end
    yStart : float, optional
        y-axis start
    yEnd : float, optional
        y-axis end
    """
    # min/max of axes
    subplot.set_xlim([xStart,xEnd])
    if not None in (yStart, yEnd):
        subplot.set_ylim([yStart,yEnd])

    # ticks
    subplot.tick_params(direction='out')
    subplot.spines['top'].set_visible(False)
    subplot.spines['right'].set_visible(False)
    subplot.get_xaxis().tick_bottom()
    subplot.get_yaxis().tick_left()

    # rotate x-labels
    xlabels = subplot.get_xticklabels()
    for label in xlabels:
        label.set_rotation(90) 
    

def plot_1d_kde_histogram(kde, subplot, 
                          xStart=1.66, xEnd=1.78, xStep=0.001,
                          xlab=None, ylab=None, title=None): 
    """Create a subplot histogram of 1d-KDE.
    
    Parameters
    ----------
    kde : KDE object
    subplot : matplotlib subplot object
    """
    xi = np.arange(xStart,xEnd,xStep)

    if kde is not None:
        yi = kde(xi)
        subplot.bar(xi, yi, width=xStep)
    
    if title:
        subplot.set_title(title[:25])
    if xlab:
        subplot.xlabel(xlab)
    if ylab:
        subplot.ylabel(ylab)

    # pretty axes
    format_subplot_axes(subplot, xStart, xEnd)


def plot_2d_kde_heatmap(kde, subplot, 
                        xStart=1.66, xEnd=1.78, xStep=0.001, 
                        yStart=1, yEnd=15000, yStep=100, 
                        xlab=None, ylab=None, title=None): 
    """Create a subplot heatmap of a 2d-KDE object.

    Parameters
    ----------
    kde : KDE object
    subplot : matplotlib subplot object
    """

    i = np.arange(xStart, xEnd, xStep)
    i_nbins = len(i)
    j = np.arange(yStart, yEnd, yStep)
    j_nbins = len(j)
    
    xi,yi = np.mgrid[i.min():i.max():i_nbins*1j, j.min():j.max():j_nbins*1j]
    
    if kde is not None:
        zi = kde(np.vstack([xi.flatten(), yi.flatten()]))
        subplot.pcolormesh(xi, yi, zi.reshape(xi.shape))
        
    if title:
        subplot.set_title(title[:25])
    if xlab:
        subplot.xlabel(xlab)
    if ylab:
        subplot.ylabel(ylab)

    # pretty axes
    format_subplot_axes(subplot, xStart, xEnd, yStart, yEnd)
        

def make_subplot_coords(n_subplot, ncol):
    """Make subplot coords.

    Parameters
    ----------
    n_subplot : int
        number of subplots to make
    ncol : int
        Number of columns with subplots
    Returns
    x : tuple 
       (subplot_coords_x-y,nrow,ncol)
    """
    # setting up subplots
    ncol = int(ncol)
    n_subplot = int(n_subplot)
    nrow = int(math.ceil(n_subplot / float(ncol)))
    
    # subplot coords 
    coords = itertools.product(range(nrow), range(ncol))
    
    return (coords, nrow)


def get_KDEs_min_max(KDEs):
    """Get the min & max of each KDE object.

    Parameters
    ----------
    KDEs : dict
        dict of KDE objects
    """
    x_min = []
    y_min = []
    x_max = []
    y_max = []    

    for KDE in KDEs.values():
        if KDE is None:
            continue            

        x_min.append(min(KDE.dataset[0]))
        x_max.append(max(KDE.dataset[0]))            
        
        if KDE.d == 2:
            y_min.append(min(KDE.dataset[1]))
            y_max.append(max(KDE.dataset[1]))
        else:
            y_min = [None]
            y_max = [None]

    return (np.min(x_min), np.max(x_max), 
            np.min(y_min), np.max(y_max))
    


def make_kde_fig(KDEs, outFile, n_subplot=0, ncol=3, 
                 xMin='', xMax='', yMin='', yMax='',
                 xX=4, yX=3.5, xStep=0.0005, yStep=100, logY=''):
    """Make KDE figure.
    
    Parameters
    ----------
    KDEs : dict
        dict of KDE objects
    outFile : str
        Name of output file.
    n_subplot : int
        Number of subplots.
    ncol : int
        Number of subplot columns.
    xMin, xMax, xStep : float
        x-axis range values
    yMin, yMax, yStep : float
        y-axis range values
    xX : float
        figure width
    yX : float
        figure height
    logY : float
        Plot y-axis as log with base of ``logY``        
    """

    # value format
    n_subplot = int(n_subplot)
    ncol = int(ncol)
    xX = float(xX)
    yX = float(yX)
    xStep = float(xStep)
    yStep = float(yStep)
    
    # checking n-dims of KDEs
    ndims = KDE_ndims(KDEs)
    if len(ndims) > 1:
        msg = 'KDEs do not all have the same number of dimensions'
    else:
        try:
            ndims = list(ndims)[0]
        except IndexError:
            msg = 'Nothing to plot!'
            raise IOError, msg

    # determine x,y min/max for plotting
    x_min,x_max,y_min,y_max = get_KDEs_min_max(KDEs)    
    if not xMin == '':
        x_min = float(xMin)
    if not xMax == '':
        x_max = float(xMax)
    if not yMin == '':
        y_min = float(yMin)
    if not yMax == '':
        y_max = float(yMax)    

    # subplot coords
    if n_subplot < 1:
        n_subplot = len(KDEs.keys())
    coords,nrow = make_subplot_coords(n_subplot, ncol)    

    # making subplot object
    f, axarr = plt.subplots(nrow, ncol, figsize=(ncol*xX, nrow*yX))
    try:         
        if axarr.ndim == 1:
            axarr = np.array([axarr])                    
    except AttributeError:
        try:
            axarr = np.array([[x] for x in axarr])
        except TypeError:
            axarr = np.array([[axarr]])            
    finally:
        if ncol == 1:
            axarr = axarr.transpose()

    # making subplots
    taxon_cnt = 0
    for (i,j),(taxon) in zip(coords, sorted(KDEs.keys())):
        msg = 'Processing KDE for taxon: "{}"\n'
        sys.stderr.write(msg.format(taxon))

        kde = KDEs[taxon]        
        subplot = axarr[i,j]
        if not logY == '':
            subplot.set_yscale('log', basey=float(logY))

        if ndims == 1:
            plot_1d_kde_histogram(kde, subplot,
                                  xStart=x_min, xEnd=x_max, xStep=xStep,
                                  title=taxon)
        elif ndims == 2:
            plot_2d_kde_heatmap(kde, subplot,
                                xStart=x_min, xEnd=x_max, xStep=xStep,
                                yStart=y_min, yEnd=y_max, yStep=yStep,
                                title=taxon)
        else:
            msg = 'KDE ndims > 2; cannot create figure'
            raise ValueError, msg

        # iter count
        taxon_cnt += 1
        if taxon_cnt >= n_subplot:
            break

    # saving figure
    plt.tight_layout()
    plt.savefig(outFile, bbox_inches='tight')
    msg = 'File written: "{}"\n'
    sys.stderr.write(msg.format(outFile))
    

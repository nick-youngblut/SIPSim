#!/usr/bin/env python

#--- Option parsing ---#
"""
diffusion: incorporate gradient diffusion in G+C variance. 

Usage:
  fragments [options] <fragment_kde>
  fragments -h | --help
  fragments --version

Options:
  <fragment_kde>    Output from the fragment_kde subcommand.
                    ('-' if input from STDIN) 
  -n=<n>            Number of Monte Carlo replicates to estimate
                    G+C error due to diffusion. 
                    [default: 10000]
  -T=<T>            Ultracentrifugation run temperature (kelvin).
                    [default: 298]
  -B=<B>            Beta coefficient based on gradient salt density.
                    [default: 1.195e9] 
  -G=<G>            G coefficient (see Clay et al., 2003). 
                    [default: 7.87e-10]
  --bw=<bw>         The bandwidth scalar or function passed to
                    scipy.stats.gaussian_kde().
  --np=<np>         Number of parallel processes.
                    [default: 1]
  --cs=<cs>         Chunksize for each process (number of taxa).
                    [default: 1]
  -h --help         Show this screen.
  --version         Show version.
  --debug           Debug mode (no parallel processes)

Description:
  The location of DNA fragments in an isopycnic gradient is
  determined by the opposing forces of sedimentation and diffusion.
  
  Sedimentation during ultracentriguation will cause a 
  homogenous population of fragments to migrate to the same
  point in the gradient (the isopycnic point); however, diffusion
  opposes this sedimentation. At equilibrium, sedimentation and
  diffusion are balanced, with no net flux of particles. Particles
  at equilibrium will be dispersed in a band whose thickness
  is dependent on diffusion. 

  Diffusion can be modeled as a gaussian point spread function.
  Calculations based on:
    Clay O, Carels N, Douady CJ, Bernardi G. (2006). Density Gradient 
    Ultracentrifugationand Whole Genome Sequences:Fine-tuning the 
    Correspondence. In:Analytical Ultracentrifugation VIII, Wandrey,
    C & Colfen, H (eds) Progress in Colloid and Polymer Science, 
    Springer Berlin Heidelberg, pp. 97-107.

  The error in 'true' G+C values caused by diffusion is estimated
  by Monte Carlo simulation. 

  **Output**
  A pickled python object {taxon_name:buoyant_density_kde} is
  written to STDOUT.
"""

# import
## batteries
from docopt import docopt
import sys,os
import cPickle as pickle 
from functools import partial
## 3rd party
import parmap
import scipy.stats as stats
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from FragGC import Frag_multiKDE
import SIPSimCython as SSC
import Utils


# functions
def make_kde(taxon_name, kde, n, bw_method, **kwargs):
    if kde is None:
        return (taxon_name, None)

    kdeBD = stats.gaussian_kde( 
        SSC.add_diffusion(
            kde.resample(size=n),
            **kwargs), 
        bw_method=bw_method)
    return (taxon_name, kdeBD)
    

def main(args):    
    """
    load kde object
    foreach kde (parallel):
      foreach iter:
        draw from kde and calculate diffusion
        append to array
      calculate BD from kde
      create kde from BD array
      return kde object
    """
    try:
        args['--bw'] =float(args['--bw'])
    except TypeError:
        pass 

    kde2d = Utils.load_kde(args['<fragment_kde>'])

    pfunc = partial(make_kde, 
                    n = int(args['-n']),
                    T = float(args['-T']),
                    B = float(args['-B']),
                    G = float(args['-G']),
                    bw_method=args['--bw'])
    
    KDE_BD = parmap.starmap(pfunc, kde2d.items(),
                            processes = int(args['--np']),
                            chunksize = int(args['--cs']),
                            parallel = not args['--debug'])
    

    pickle.dump({taxon:KDE for taxon,KDE in KDE_BD}, 
                sys.stdout)
    
    
# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')            
    main(args)
    

        

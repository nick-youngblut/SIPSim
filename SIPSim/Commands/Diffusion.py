#!/usr/bin/env python

"""
diffusion: incorporate gradient diffusion error in fragment buoyant density estimates

Usage:
  diffusion [options] <fragment_kde>
  diffusion -h | --help
  fragments --version

Options:
  <fragment_kde>      Output from the fragment_kde subcommand.
                      ('-' if input from STDIN) 
  -m=<m>              Method for calculating diffusion error
                      ('Clay' or 'Meselson').
                      [Default: Clay]
  -B=<B>              Beta coefficient based on gradient salt density.
                      [default: 1.14e9]
  -D=<D>              Average particle density in gradient.
                      [default: 1.7]
  -w=<w>              Angular velocity of rotor (omega^2).
                      [default: 33172837]
  --r_min=<rm>        radius min from axis of rotation (cm).
                      [default: 2.6]
  --r_max=<rx>        radius max from axis of rotation (cm).
                      [default: 4.85]
  -t=<t>              Ultracentrifugation run time (seconds).
                      [default: 23760]
  -T=<T>              Ultracentrifugation run temperature (kelvin).
                      [default: 293.15]
  -G=<G>              G coefficient (see Clay et al., 2003). 
                      [default: 7.87e-10]
  -M=<M>              Molecular weight per base pair of dry cesium DNA.
                      [default: 882]
  -n=<n>              Number of Monte Carlo replicates to estimate
                      G+C error due to diffusion. 
                      [default: 500000]
  --BD_range=<bdr>    Start,stop,step for making fragment BD bins.
                      [default: 1.66,1.78,0.001]
  --len_range=<lr>    Start,stop,step for making fragment length bins.
                      [default: 200,75000,100]
  --bw=<bw>           The bandwidth scalar or function passed to
                      scipy.stats.gaussian_kde().
  --np=<np>           Number of parallel processes.
                      [default: 1]
  -h --help           Show this screen.
  --version           Show version.
  --debug             Debug mode (no parallel processes)

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
    Clay O, Douady CJ, Carels N, Hughes S, Bucciarelli G, Bernardi G. (2003). 
    Using analytical ultracentrifugation to study compositional variation in 
    vertebrate genomes. Eur Biophys J 32: 418-426.

    Meselson M, Stahl FW, Vinograd J. (1957). Equilibrium Sedimentation
    of Macromolecules in Density Gradients. PNAS 43: 581-588.

  The error in 'true' G+C values caused by diffusion is estimated
  by Monte Carlo simulation. 

  Output
  ------
  A pickled python object {taxon_name:buoyant_density_kde} is
  written to STDOUT.
"""

# import
## batteries
from docopt import docopt
import sys,os
## application libraries
from SIPSim import Diffusion

        
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    args['-n'] = float(args['-n'])
    Diffusion.main(args)
   

        

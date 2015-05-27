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
                    [default: 100000]
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
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import Diffusion

    
    
# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')            
    Diffusion.main(args)
    

        

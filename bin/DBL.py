#!/usr/bin/env python

#--- Option parsing ---#
"""
DBL: include diffusion boundary layer contamination into G+C variance

Usage:
  DBL [options] <fragment_kde>
  DBL -h | --help
  DBL --version

Options:
  <fragment_kde>      Output from the fragment_kde subcommand.
                      ('-' if input from STDIN) 
  -B=<B>              Beta coefficient.
                      [default: 1.14e9]
  -D=<D>              Average particle density in gradient.
                      [default: 1.7]
  -w=<w>              Angular velocity of rotor (omega^2).
                      [default: 33172781]
  --tube_diam=<td>    cfg tube diameter (cm).
                      [default: 1.3]
  --tube_height=<th>  cfg tube height (cm).
                      [default: 4.8]
  --r_min=<rm>        radius min from axis of rotation (cm).
                      [default: 2.6]
  --r_max=<rx>        radius max from axis of rotation (cm).
                      [default: 4.85]
  --vertical          Vertical rotor (instead of fixed angle).
  --frac_abs=<fa>     Fraction of DNA absorbed to the cfg tube wall.
                      [default: 0.001]
  --BD_min=<bm>       Min BD used to determine the DBL.
                      [default: 1.59]
  --BD_max=<bx>       Max BD used to determine the DBL.
                      [default: 1.77]
  --DBL_out=<z>       Write the DBL index to file named <z>.
                      [default: DBL_index.txt]
  -n=<n>              Number of Monte Carlo replicates to estimate
                      G+C error due to DBL. 
                      [default: 500000]
  --comm=<c>          Community file. Used for scaling fraction absorbed
                      by the pre-fraction abundance.
                      Fraction absorbed will be constant if `comm` not provieded.
  --commx=<cx>        Scaling factor for fraction_absorbed ~ pre-fraction_abund.
                      The intercept is set by the `frac_abs` option.
                      [default: 0.8]
  --bw=<bw>           The bandwidth scalar or function passed to
                      scipy.stats.gaussian_kde(). 
  --np=<np>           Number of parallel processes.
                      [default: 1]
  -h --help           Show this screen.
  --version           Show version.
  --debug             Debug mode (no parallel processes)

Description:
  In isopycnic centrifugation, the gradient forms perpendicular to the
  axis of rotation. Thus, in a fixed angle or vertical rotor, the gradient 
  is oriented no down the length of the cfg tube, but at an angle. 
  The gradient formed during ultracentrifugation will rotate once the tube
  is placed in a vertical orientation for fractionation.
  For more info, see: 
      Flamm WG, Bond HE, Burr HE. (1966). 
      Density-Gradient centrifugation of DNA in a fixed-angle rotor. 
      Biochimica et Biophysica Acta (BBA) -
      Nucleic Acids and Protein Synthesis 129: 310-317.

  This difference in gradient orientation between centrifugation and 
  fractionation may result in gradient 'smearing' if DNA binds to the tube
  wall during centrifugation and then diffuses back into solution during the 
  period between centrigution and fractionation. This diffusive boundary
  layer (DBL) can introduce 'light' DNA into 'heavy' fractions and vice versa.
  For instance:
    In a fixed angle rotor, particles with a buoyant density of 1.7 (at
    equilibrium) will be collected in a band spanning a couple of centimeters
    of the tube height.
           /  /
    top-->/  /
         /| /     <-- cfg tube during ultracentrifugation 
        / |/
       (  /<--radient band bottom (BD of 1.7)
        ^^

      |  |
      |  |
      |--| <--gradient band (BD of 1.7) during fractionation
      \__/

  Based on information on the rotor and cfg run conditions, this script
  models the gradient 'smearing' that would result from a diffusive boundary
  layer (contaminating DNA diffusing back into the gradient from the tube wall).

  The error in 'true' G+C values caused by the DBL is estimated
  by Monte Carlo simulation. 

  The default rotor & tube parameters are for a Beckman tubes (ref# 361621) 
  in a Beckman TLA-110 rotor spinning at 55000 rpm. 

  Vertical rotors: use the `vertical` option.

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
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

import DiffBoundLayer

        
# main
if __name__ == '__main__':
    args = docopt(__doc__, version='0.1')            
    args['-n'] = float(args['-n'])
    DiffBoundLayer.main(args)
    

        

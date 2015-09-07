#!/usr/bin/env python

#--- Option parsing ---#
"""
PCR: simulate PCR of gradient fraction DNA samples

Usage:
  PCR [options] <OTU_table>
  PCR -h | --help
  PCR --version

Options:
  <OTU_table>              OTU table file.
  --n_cycles=<n>           Number of PCR cycles.
                           [Default: 30]
  --DNA_conc_dist=<dc>     Distribution of starting DNA molarities (units = uM).
                           Use 'uniform' if all reactions used the same amount
                           of input DNA.
                           (see numpy.random for a list of distributions)
                           [Default: uniform]
  --DNA_conc_dist_p=<dp>   Distribution parameters.
                           (see numpy.random for a list of parameters)
                           [Default: low:0.3,high:0.3]
  --primer_conc=<pc>       Molarity of forward and reverse primers (units = uM).
                           [Default: 1]
  -f=<f>                   The initial PCR reaction effiency.
                           [Default: 1]
  -k=<k>                   k parameter used in Suzuki & Giovannoni (1996).
                           [Default: 5]
  --version                Show version.
  --debug                  Debug mode (turns off parallel processing)
  -h --help                Show this screen.


Description:
  Simulate PCR on the template DNA for each gradient fraction sample.
  This simulation will account for template saturation, where 
  PCR effeciency declines with increased template concentrations in later
  PCR cycles (see Suzuki & Giovannoni, 1996).

References:
  Suzuki MT, Giovannoni SJ. (1996). Bias caused by template annealing in the
  amplification of mixtures of 16S rRNA genes by PCR. Appl Environ Microbiol
  62:625-630.
"""

# import
## batteries
from docopt import docopt
import os, sys
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)

from Utils import parseKeyValueString as distParamParse
from OTU_Table import OTU_table
from PCR_Sim import PCR_Sim
    

def main(Uargs):
    # parsing dist params
    Uargs['--DNA_conc_dist_p'] = distParamParse(Uargs['--DNA_conc_dist_p'])

    # loading OTU table 
    otu_tbl = OTU_table.from_csv(Uargs['<OTU_table>'], sep='\t')

    # PCR simulation
    otu_tbl_pcr = PCR_Sim(otu_tbl,
                          DNA_conc_dist = Uargs['--DNA_conc_dist'],
                          DNA_conc_dist_p = Uargs['--DNA_conc_dist_p'],
                          primer_conc = float(Uargs['--primer_conc']),
                          n_cycles = int(Uargs['--n_cycles']),
                          f_0 = float(Uargs['-f']),
                          k = float(Uargs['-k']))
    

# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    main(Uargs)

    

        

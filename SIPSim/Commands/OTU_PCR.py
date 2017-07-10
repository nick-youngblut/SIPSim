#!/usr/bin/env python

"""
OTU_PCR: simulate PCR of gradient fraction DNA samples

Usage:
  OTU_PCR [options] <OTU_table>
  OTU_PCR -h | --help
  OTU_PCR --version

Options:
  <OTU_table>              OTU table file.
  --n_cycles=<n>           Number of PCR cycles.
                           [Default: 30]
  --DNA_conc_dist=<dc>     Distribution of starting DNA molarities for each
                           sample (units = uM).
                           Note: Use 'uniform' if all reactions used the same
                           amount of input DNA.
                           [Default: uniform]
  --DNA_conc_dist_p=<dp>   Distribution parameters.
                           (see numpy.random for a list of parameters)
                           [Default: low:0.3,high:0.3]
  --primer_conc=<pc>       Molarity of forward and reverse primers (units = uM).
                           [Default: 1]
  --ratio=<r>              Amplicon to primer length ratio.
                           [Default: 10]
  -f=<f>                   The theoretical maximum PCR efficiency.
                           [Default: 1]
  -k=<k>                   k parameter used in Suzuki & Giovannoni (1996).
                           [Default: 5]
  --version                Show version.
  --debug                  Debug mode (verbose output)
  -h --help                Show this screen.


Description:
  Simulate PCR on the template DNA for each gradient fraction sample.

  This simulation will account for template saturation, where 
  PCR effeciency declines with increased template concentrations in later
  PCR cycles (see Suzuki & Giovannoni, 1996).

  Output
  ------
  A tab-delimited OTU table written to STDOUT.

References:
  Suzuki MT, Giovannoni SJ. (1996). Bias caused by template annealing in the
  amplification of mixtures of 16S rRNA genes by PCR. Appl Environ Microbiol
  62:625-630.
"""

# import
## batteries
from docopt import docopt
import sys
import os
## application libraries
from SIPSim.Utils import parseKeyValueString as distParamParse
from SIPSim.OTU_Table import OTU_table
from SIPSim.PCR import PCR_sim
    

def main(args):
    # parsing dist params
    args['--DNA_conc_dist_p'] = distParamParse(args['--DNA_conc_dist_p'])

    # loading OTU table 
    otu_tbl = OTU_table.from_csv(args['<OTU_table>'], sep='\t')

    # PCR simulation
    PCR_sim(otu_tbl,
            DNA_conc_dist = args['--DNA_conc_dist'],
            DNA_conc_dist_p = args['--DNA_conc_dist_p'],
            primer_conc = float(args['--primer_conc']),
            n_cycles = int(args['--n_cycles']),
            f_0 = float(args['-f']),
            k = float(args['-k']),
            ratio = float(args['--ratio']),
            debug=args['--debug'])

    # writing out file
    otu_tbl.sort_values(by=['library','taxon','BD_mid'], inplace=True)
    otu_tbl.to_csv(sys.stdout, sep='\t', index=False)
    


def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)

    

        

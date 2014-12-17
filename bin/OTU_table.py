#!/usr/bin/env python

#--- Option parsing ---#
"""
OTU_table: create an OTU table of gradient fractions based on simulated
fragment G+C content and isotope incorporation

Usage:
  OTU_table [options] <fragGC_file> <comm_file> <incorp_file> <frac_file>
  OTU_table -h | --help
  OTU_table --version

Options:
  <fragGC_file>       Name of file produced by fragGC subcommand.
  <comm_file>         Name of file produced by gradientComms subcommand.
  <incorp_file>       Name of file produced by isoIncorp subcommand.
  <frac_file>         Name of file produced by fractions subcommand.
  --abs_abund=<aa>    Absolute abundance of all taxa in the community. [default: 1e6]
  --g_noise=<gn>      scipy distribution function describing gradient 'noise'. [default: cauchy]
  --scale=<np>        Scale parameter for noise gradient distribution. [default: 0.0001]
  --a_weight=<aw>     Abundance weighting for isotope incorporation.
  --isotope=<is>      Isotope incorporated by taxa (13C or 15N). [default: 13C]
  -h --help           Show this screen.
  --version           Show version.
  --debug             Debug mode

Description:

"""

# import
## batteries
from docopt import docopt
import os, sys
import logging
from collections import defaultdict
import time
## 3rd party
import pandas as pd
## application libraries
scriptDir = os.path.dirname(__file__)
libDir = os.path.join(scriptDir, '../lib/')
sys.path.append(libDir)


from SIPSim import FragGC_KDE, CommTable, FracTable, IsoIncorpTable, OTU_table 


# functinos
def main(Uargs):
    
    # --abs_abund as int
    try:
        Uargs['--abs_abund'] = int(float(Uargs['--abs_abund']))
    except KeyError:
        raise KeyError('Cannot find "--abs_abund" key')
    except TypeError:
        raise TypeError('"{}" must be float-like'.format(str(Uargs['--abs_abund'])))
    
    
    # loading fragGC file
    fragGC_kde = FragGC_KDE(Uargs['<fragGC_file>'])

    # loading community file
    comm = CommTable.from_csv(Uargs['<comm_file>'], sep='\t')
    comm.set_abs_abund(Uargs['--abs_abund'])
    
    
    # loading incorp file
    incorp = IsoIncorpTable.from_csv(Uargs['<incorp_file>'], sep='\t')

    # debug
    #for (x,y,z) in incorp.iter_incorpFuncs():
    #    print '{}->{}:\n{}'.format(x,y,z)
    #sys.exit();
    
    
    # loading fraction file
    frac = FracTable.from_csv(Uargs['<frac_file>'], sep='\t')

    
    # initializing OTU table class
    OTU = OTU_table(frac,
                    Uargs['--g_noise'],
                    Uargs['--scale'],
                    Uargs['--a_weight'],                    
                    Uargs['--isotope'])

    
    # checking on library overlap
    if not OTU.checkLibOverlap([
            [x for x in comm.iter_libraries()],
            [x for x in incorp.iter_libraries()],
            [x for x in frac.iter_libraries()]]):
        logging.warning('Not all tables contain the same library IDs')

        
    # iter by library:
    isotopeMaxBD = OTU.get_isotopeMaxBD()
    OTU_counts = dict()
    for libID in comm.iter_libraries():
        sys.stderr.write('Processing library: "{}"\n'.format(libID))
        
        # make dataframe for OTU counts: taxa X fractions
        OTU_counts[libID] = OTU.make_emptyCountTable(comm.get_unique_taxon_names(),
                                                     frac.get_libFracIDs(libID))

        # all values to integers
        OTU_counts[libID] = OTU_counts[libID].astype(int)

        # iter by taxon:
        for taxon_name in comm.iter_taxa(libID):
            sys.stderr.write('  Processing taxon: "{}"\n'.format(taxon_name))
            
            taxonAbsAbund = comm.get_taxonAbund(libID, taxon_name)

            # sampling GC value from taxon KDE
            t0 = time.time()
            GC_vals = fragGC_kde.sampleTaxonKDE(taxon_name, n_samples=taxonAbsAbund)

            # sampling intra-taxon incorp for taxon; return: iterator
            t1 = time.time()
            incorp_vals = incorp.sample_incorpFunc(libID, taxon_name, n_samples=taxonAbsAbund)
            
            # iter GC value:
            t2 = time.time()
            for gc_val in GC_vals:
                # raw BD based on GC
                BD = gc_val / 100 * 0.098 + 1.66

                # BD + BD shift from isotope incorporation
                ## TODO: implement abundance-weighting
                incorp_perc = incorp_vals.next()
                BD = BD + isotopeMaxBD * (incorp_perc / 100)

                # simulating gradient noise
                BD = OTU.sample_g_noise_func(BD)[0]
                
                # determine the fraction that would contain the fragment
                fracID = frac.which_frac(libID, BD)

                # adding to OTU count table
                try:
                    OTU_counts[libID][fracID].loc[taxon_name] += 1
                except ValueError:
                    pass

            t3 = time.time()
            #print [t1 - t0, t2 - t1, t3 - t2]; sys.exit()
                    
        # replace fraction IDs in column names with BD_min-BD_max
        OTU_counts[libID].columns = frac.fracID2BDminmax(libID, OTU_counts[libID].columns)

                    
    # combine the dataframes and write
    pd.concat(OTU_counts.values(), axis=1).to_csv(sys.stdout, sep='\t')


        
# main
if __name__ == '__main__':
    Uargs = docopt(__doc__, version='0.1')
    main(Uargs)
    

        
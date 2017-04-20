import sys
import os
import pandas as pd
import numpy as np
import scipy


class fraction(object):
    """class for simulating fractions"""
    
    def __init__(self):
        pass
#        self.nGradient = nGradient

    def simFractions(self, BD_min=1.66, BD_max=1.76, fracSizeMean=0.003, fracSizeStdev=0.001, fracSizeMin=0.001):
        """Simulation gradient fractions. Drawing from a normal distribution

        Args:
        BD_min -- min BD in gradient
        BD_max -- max BD in gradient
        fracSizeMean -- mean fraction size
        fracSizeStdev -- stdev of fraction sizes
        """

        # simulating fraction sizes
        fracRange = BD_max - BD_min
        fracs = []
        fracSum = 0
        while 1:
            fracSize = np.random.normal(loc=fracSizeMean, scale=fracSizeStdev)
            if fracSize <= 0 or fracSize < fracSizeMin:
                continue
            elif fracSum + fracSize > fracRange:
                continue
            else:
                fracs.append(fracSize)
                fracSum = sum(fracs)

            # end loop?
            if fracSum >= fracRange:
                break
            #elif fracSum > fracRange:
            #    sys.exit('Internal Error')
            else:
                continue

        # making fraction start-end DataFrame
        tbl = dict(start=[], end=[])
        [tbl['start'].append(x + BD_min) for x in fracs]
        [tbl['end'].append(x + BD_min) for x in fracs]
        
        self.fracs = pd.DataFrame(tbl)
        
                                    
        
    def simGradientFracs(self, nGradients=1):
        pass


#-- unit tests --#
if __name__ == '__main__':
    frac = fraction()    
    print frac.simFractions()
    
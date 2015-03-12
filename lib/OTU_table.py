import os, sys
import pandas as pd
import numpy as np

from Utils import _table


class OTU_table(_table):

    def __init__(self, *args, **kwargs):
        _table.__init__(self, *args, **kwargs)

        
#    def wide2long(self, id_vars=['library', 'fractions'],
#                  var_name='taxon', value_name='count'):
#        """Convert table from wide to long format.
#        """
#        self.df = pd.melt(self.df,
#                id_vars=id_vars,
#                var_name=var_name,
#                value_name=value_name)        
#        print self.df;
#        sys.exit()
        

#    def long2wide(self):
#        """Convert table from long to wide format.
#        """
#        self.df = pd.pivot_table(self.df, values='count',
#                       index='taxon', columns=['library','fractions'])
#        self.df.columns =  ['__'.join(x) for x in self.df.columns.tolist()]

        
    def write_table(self, sep='\t', outBuf=sys.stdout, index=False):
        """Writing out table.
        Args:
        sep -- column separator
        outBuf -- output buffer object
        """
        self.df.to_csv(outBuf, sep=sep, index=index)
        
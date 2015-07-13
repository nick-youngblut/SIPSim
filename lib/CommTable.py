from Utils import _table
import sys, os
import re
import numpy as np

class CommTable(_table):
    """pandas DataFrame subclass
    """
    
    def __init__(self, *args, **kwargs):
        """Subclassing pandas dataframe. CommTable.from_[csv,etc]()
        Can be used for reading in a table.
        """
        _table.__init__(self, *args, **kwargs)


    def taxonInLib(self, taxon_name, libID):
        """Checking whether taxon in the selected library
        Args:
        taxon_name -- taxon name
        libID -- library ID
        Return:
        boolean
        """
        libID = re.sub('\D+', '', libID)
        dfSub = self.df.loc[self.df['taxon_name'] == taxon_name]
        libs = [str(x) for x in dfSub['library'].tolist()]
        return libID in libs

            
    def get_taxonAbund(self, taxon_name, libID=None, abs_abund=False):
        """Getting the abundance(s) of taxon in the comm file.
        Args:
        taxon_name -- name of taxon
        libID -- library ID. If None, all libraries selected.
        abs_abund -- return absolute abundance instead of relative abundance
        Return:
        iterable -- [abundance values for the taxon]
        """
        retCol = 'abs_abund' if abs_abund else 'rel_abund_perc'
        assert retCol in self.df.columns, \
            '"{}" column not found'.format(retCol)
        if libID is not None:            
            df_sub =  self.df.loc[(self.df['library'] == libID) &
                                  (self.df['taxon_name'] == taxon_name)]
        else:
            df_sub =  self.df.loc[(self.df['taxon_name'] == taxon_name)]

        return df_sub[retCol].tolist()

    
    def get_unique_libIDs(self):
        """Getting all unique libIDs from the community table
        """
        return self.df['library'].unique().tolist()

    def get_unique_taxon_names(self):
        """Getting all unique taxon names from the community table
        """
        return self.df['taxon_name'].unique()        

    @property
    def abs_abund(self):
        return self.df['abs_abund']

    @abs_abund.setter
    def abs_abund(self, abs_abund):
        """Setting the absolute abundance of each taxon based on user-value.
        Args:
        abs_abund -- int; total abundance of all taxa in each library
        Returns:
        None
        """
        try:
            x = self.df['rel_abund_perc']
            self.df['abs_abund'] = x / 100 * int(abs_abund)
        except KeyError:
            raise KeyError('"rel_abund_perc" column not found in comm file')

        self.df['abs_abund'] = np.round(self.df['abs_abund'], 0)

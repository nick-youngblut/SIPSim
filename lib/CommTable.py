from Utils import _table
import re

class CommTable(_table):
    
    def __init__(self, *args, **kwargs):
        """Subclassing pandas dataframe. CommTable.from_[csv,etc]()
        Can be used for reading in a table.
        """
        _table.__init__(self, *args, **kwargs)

        
    def set_abs_abund(self, abs_abund):
        """Setting the absolute abundance of each taxon based on user-value.
        Args:
        abs_abund -- int; total abundance of all taxa in each library
        """
        try:
            self.df['abs_abund'] = self.df['rel_abund_perc'] / 100 * int(abs_abund)
        except KeyError:
            raise KeyError('"rel_abund_perc" column not found in comm file')


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

            
    def get_taxonAbund(self, taxon_name, libID=None):
        """Getting the abundance(s) of taxon in the comm file.
        Args:
        taxon_name -- name of taxon
        libID -- library ID. If None, all libraries selected        
        Return:
        iterable of abundance values for the taxon
        """
        retCol = 'rel_abund_perc'                    
        assert retCol in self.df.columns, \
            '"{}" column not found'.format(retCol)
        if libID is not None:            
            df_sub =  self.df.loc[(self.df['library'] == libID) &
                                  (self.df['taxon_name'] == taxon_name)]
        else:
            df_sub =  self.df.loc[(self.df['taxon_name'] == taxon_name)]
        
        return df_sub[retCol].tolist()

    
    def get_unique_libIDs(self):
        return self.df['library'].unique().tolist()

    def get_unique_taxon_names(self):
        return self.df['taxon_name'].unique()        

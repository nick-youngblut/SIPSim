from Utils import _table

class CommTable(_table):
    
    def __init__(self, *args, **kwargs):
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

            
    def get_taxonAbund(self, libID, taxon_name, rel=False):
        """Getting the abundance of taxon in the comm file.
        Absolute abundance is returned unless rel=True.
        If rel=True, relative abundance is returned.
        Args:
        rel -- return relative abundance instead of absolute?
        """
        if rel == True:
            retCol = 'rel_abund_perc'
        else:
            retCol = 'abs_abund'
            
        assert retCol in self.df.columns, \
            '"{}" column not found'.format(retCol)
        df_sub =  self.df.loc[(self.df['library'] == libID) &
                              (self.df['taxon_name'] == taxon_name)]
        return int(df_sub[retCol])

    
    def get_unique_taxon_names(self):
        return self.df['taxon_name'].unique()        

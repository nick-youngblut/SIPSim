# import
## batteries
import sys
## 3rd party
from intervaltree import Interval, IntervalTree
## application
from Utils import _table
from OrderedSet import OrderedSet
        

class FracTable(_table):
    """Subclass of pandas.DataFrame

    TODO
    ----
    Assertion that fractions are ordered and not overlapping.
    """
    
    def __init__(self, *args, **kwargs):
        """_table subclass"""
        _table.__init__(self, *args, **kwargs)
        
        # 'delete' inherited methods involving taxa
        iter_taxa = property()
        iter_taxonRowsInLib = property()

        # assertions on columns
        assertCols = ['library','fraction','BD_min','BD_max','fraction_size']
        for col in assertCols:
            msg = 'Column "{}" not found in file "{}"'
            assert col in self.df.columns, msg.format(col, self.tableFileName)

        # fraction ID as int
        self.df['fraction'] = self.df['fraction'].astype(int)

        # creating an interval tree of fractions
        self._set_itrees()

        
    def _set_itrees(self):
        """Setting a dict of interval trees for the BD-min/max ranges
        (1 per lib). Max in interval ranges is non-inclusive.
        in-place edit of ``trees`` attribute.
        """
        self.itrees = dict()

        def make_frac_interval(x):
            # range start, range end, fraction ID
            return Interval(float(x[0]),float(x[1]),int(x[2]))
        
        for libID in self.iter_libraries():
            x = self.df['library'] == libID, ['BD_min','BD_max','fraction']
            df_sub = self.df.loc[x]
            y = list(df_sub.apply(make_frac_interval, axis=1))
            self.itrees[libID] = IntervalTree(y)

            
    def which_frac(self, libID, BD_value):
        """Determine which of the simulated fractions that the BD value falls
        into. Using an interval tree of BD-min/max values.
        BD_max is non-inclusive.
        
        Parameters
        ----------
        libID : str
            library ID
        BD_value : float
            Bouyant density value

        Returns
        -------
        fracIDs : list 
            fractionIDs that the BD_value falls into
        """
        libID = str(libID)
        BD_value = float(BD_value)
        
        assert hasattr(self, 'itrees'), 'Cannot find itrees attrib'
        
        try:
            fracIDs = [x.data for x in self.itrees[libID][BD_value]]            
        except KeyError:
            msg = 'Library "{}" not found in fraction table'
            raise KeyError(msg.format(libID))

        fracIDs_len = len(fracIDs)
        if fracIDs_len == 0:
            return None
        elif fracIDs_len == 1:
            return fracIDs[0]
        else:
            msg = 'BD value "{}" in library "{}" falls into >1 fraction!'
            raise ValueError(msg.format(str(BD_value), libID))
        
            
    def fracID2BDminmax(self, libID, fracIDs=None):
        """Convert fractionIDs to 'BDmin-BDmax'.
        
        Parameters
        ----------
        libID : str
            library ID
        fracIDs : iterable
            fraction IDs; if None: all fracIDs used

        Returns
        -------
        BD_min_max : list 
            A list of 'BDmin-BDmax' (1 string per fraction).
        """
        df_sub = self.df.loc[(self.df['library'] == libID)]
        df_sub.index = df_sub['fraction']

        if fracIDs is not None:
            df_sub = df_sub.reindex(fracIDs)
        else:
            pass
            
        
        def catCol(x, sep='-', cols=['BD_min','BD_max']):
            return sep.join([libID] + [str(y) for y in x[cols]])
        
        return  list(df_sub.apply(catCol, axis=1))

        
    def get_fracBDminxax(self, libID, fracID):
        """Get the BD min/max values for a fraction in a library.
        
        Parameters
        ----------
        libID : str
             library ID
        fracID :str
             fraction ID

        Returns
        -------
        BD_min_max : list
            [BD_min, BD_max]
        """
        x = (self.df['library'] == libID) & (self.df['fraction'] == fracID) 
        df_sub = self.df.loc[x]
        return list(df_sub.loc[0][['BD_min','BD_max']])


    def BD_bins(self, libID=None):
        """Return an ordered set of BD min-max bins.
        Assumptions:
        * bins do not overlap
        * bins cover the whole BD range (min[n] == max[n-1])
        Useful for binning with np.digitize().

        Parameters
        ----------
        libID : str
            library ID; select BD values from just 1 library

        Returns
        -------
        BD_bins : OrderedSet 
            BD min-max bins
        """
        if libID is not None:
            df_sub = self.df.loc[self.df['library'] == libID]
        else:
            df_sub = self.df

        x = zip(df_sub['BD_min'].tolist(), df_sub['BD_max'].tolist())
        return OrderedSet([j for i in x for j in i])

        
    def get_fracIDs(self):
        """Return list of library-fractionIDs.
        Format: [libraryID]_[fractionID]_[BD_min]_[BD_max]
        """
        df_sub = self.df[['library','fraction','BD_min','BD_max']]
        f = lambda x : '_'.join([str(y) for y in x])
        return list(df_sub.apply(f, axis=1))

    
    def get_libFracIDs(self, libID):
        """Return list of fractionIDs in a library
        
        Parameters
        ----------
        libID : str
            library ID

        Returns
        -------
        fracIDs : list 
            fraction IDs
        """
        df_sub = self.df.loc[self.df['library'] == libID]
        return list(df_sub['fraction'])
        
        
    def get_libFracDict(self):
        """Get a dict of info for each fraction in each library

        Returns
        -------
        d : dict
             {libraryID : {fractionID : dict_of_info}}            
        """
        libFrac = defaultdict(defaultdict(dict))
        for i in xrange(self.df.shape[0]):
            row = self.df.loc[i]
            libFrac[str(row['library'])][str(row['fraction'])] = dict()
        return libFrac

# import
import sys
import numpy as np
import dendropy 
from dendropy.utility.error import DataParseError


def read_tree(tree_obj):
    """Parsing a tree file that's in various possible formats
    Args:
    tree_obj -- path or data stream fo tree file
    Returns:
    dendropy object 
    """
    if tree_obj is None:
        return None

    psblFormats = ['newick', 'nexus']
    for f in psblFormats:
        try:
            tree = dendropy.Tree.get_from_path(tree_obj, schema=f)
        except (DataParseError, TypeError):
            pass
        else:
            break
        try:
            tree = dendropy.Tree.get_from_stream(tree_obj, schema=f)
        except (DataParseError, ValueError):
            pass
        else:
            break
        
        
    if isinstance(tree, dendropy.dataobject.tree.Tree):
        return tree
    else:
        msg = 'Do not recongnize tree file format (possible formats:{})'
        fmts = ','.join(psblFormats)
        raise IOError, msg.format(fmts)


def get_leaf_node_labels(tree):
    """Return a list of leaf labels.
    Args:
    tree -- dendropy tree object
    Returns:
    list -- [leaf_label, ...]
    """
    return [x.taxon for x in tree.leaf_iter()]



class BrownianMotion:
    
    def __init__(self, tree, start=0, sigma=0.1, ratio=0.5, minVal=0, 
                 maxVal=100, verbose=False):
        """Calls sim_traits function and stores tree with simulated characters
        as an attribute.
        Args:
        tree -- dendropy tree object. It will be deepcopied
        start -- starting value for continuous character evolution
        sigma -- sigma use for drawing from a normal distribution
        ratio -- weight parameter for random vs Brownian motion
        range: 0-1; 0 = purely random; 1 = purely Brownian
        minVal -- minimum bounds on simulated values
        maxVal -- maximum bounds on simulated values
        verbose -- verbose output        
        """
        if tree is None:
            msg = 'tree object cannot be None'
            raise TypeError, msg
        self.tree = dendropy.deepcopy(tree)
        self.start = start
        self.sigma = sigma
        self.ratio = ratio
        self.minVal = minVal
        self.maxVal = maxVal
        self.verbose = verbose

        self.sim_traits()


    def __repr__(self):
        s = ''
        for node in self.tree.leaf_iter():
            s += '{} => {}\n'.format(node.taxon, node.value)
        return s


    def sim_traits(self):
        """Trait simulation as detailed in:
        author = {Munkemuller, Tamara and Lavergne, Sebastien and Bzeznik,
        Bruno and Dray, Stephane and Jombart,
        Thibaut and Schiffers, Katja and Thuiller, Wilfried},
        title = {How to measure and test phylogenetic signal},
        journal = {Methods in Ecology and Evolution}        
        Returns:
        in-place edit of self.tree
        """        
        ntaxa = len(self.tree.nodes())
        # simulate brownian motion
        BM = np.random.normal(loc=0, scale=self.sigma, size=ntaxa)
        BM = np.cumsum(BM) #+ self.start
        # random values (random ordering of cumsum(BM) values)
        rnd = np.random.permutation(BM)
        # making weighted sums
        ws = self.ratio * BM + (1-self.ratio) * rnd
        # z-scaling weighted sums
        ws = (ws - np.mean(ws)) / np.std(ws)
        # adding start
        ws += self.start
 
        # setting values within bounds
        if self.minVal is not None:
            ws = [x if x >= self.minVal else self.minVal for x in ws]
        if self.maxVal is not None:
            ws = [x if x <= self.maxVal else self.maxVal for x in ws]
    
        # setting node values
        for i, node in enumerate(self.tree.preorder_node_iter()):
            node.value = ws[i]
            if self.verbose and node.taxon is not None:
                print('{} : {}'.format(node.taxon, node.value))

    
    def sample(self, taxon_label):
        """Getting value for taxon in self.tree if taxon is present;
        otherwise, returns value from random uniform with range: (0,100).
        Args:
        taxon_label -- string of a taxon in self.tree
        Returns:
        list -- [float]
        """
        func = lambda taxon: True if taxon.label == taxon_label else False
        node = self.tree.find_node_with_taxon(func)
        if node is None:
            msg = 'Taxon: "{}" not found in the provided phylogeny. Setting random value.\n'
            sys.stderr.write(msg.format(taxon_label))
            return np.random.uniform(0,100)
        else:
            return [node.value]

    

if __name__ == '__main__':
    mle = dendropy.treesim.birth_death(birth_rate=1, death_rate=0.5, ntax=10)
    sim_traits(mle, verbose=True)
    mle.print_plot(display_width=70)

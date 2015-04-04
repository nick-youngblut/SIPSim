# import
import numpy as np
import dendropy 
from dendropy.utility.error import DataParseError


def read_tree(tree_obj):
    """Parsing a tree file that's in various possible formats
    Args:
    tree_obj -- path or data stream fo tree file
    Return:
    dendropy object if successful; else raise DataParseError
    """
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
        raise DataParseError, msg.format(fmts)



def sim_traits(tree, start=0, sigma=0.1, weight=0.5, verbose=False):
    """Trait simulation as detailed in:
    author = {Munkemuller, Tamara and Lavergne, Sebastien and Bzeznik,
              Bruno and Dray, Stephane and Jombart,
              Thibaut and Schiffers, Katja and Thuiller, Wilfried},
    title = {How to measure and test phylogenetic signal},
    journal = {Methods in Ecology and Evolution}
    
    Args:
    tree -- dendropy tree object
    start -- starting value for continuous character evolution
    sigma -- sigma use for drawing from a normal distribution
    weight -- weight parameter for random vs Brownian motion
              range: 0-1; 0 = purely random; 1 = purely Brownian
    verbose -- verbose output
    """
    ntaxa = len(tree.nodes())
    # simulate brownian motion
    BM = np.random.normal(loc=0, scale=sigma, size=ntaxa)
    BM = np.cumsum(BM) + start
    # random values
    rnd = np.random.permutation(BM)
    # making weighted sums
    ws = weight * BM + (1-weight) * rnd
    # z-scaling weighted sums
    ws = (ws - np.mean(ws)) / np.std(ws)
    
    for i, node in enumerate(tree.preorder_node_iter()):
        node.value = ws[i]
        if verbose and node.taxon is not None:
            print('{} : {}'.format(node.taxon, node.value))
            


class BrownianMotion:
    
    def __init__(self, tree, start=0, sigma=0.1, weight=0.5, verbose=False):
        """Calls sim_traits function and stores tree with simulated characters
        as an attribute.
        Args:
        tree -- dendropy tree object or file with tree
        start -- starting value for continuous character evolution
        sigma -- sigma use for drawing from a normal distribution
        weight -- weight parameter for random vs Brownian motion
        range: 0-1; 0 = purely random; 1 = purely Brownian
        verbose -- verbose output        
        """
        if tree is None:
            self.tree = None
        else:
            self.tree = read_tree(tree)
        sim_traits(self.tree, start=start, sigma=sigma, 
                   weight=weight,verbose=verbose)

        
    def sample(self, taxon_label):
        assert self.tree is not None, 'No tree object to sample'
        func = lambda taxon: True if taxon.label == taxon_label else False
        node = self.tree.find_node_with_taxon(func)
        print node.value; sys.exit()


if __name__ == '__main__':
    mle = dendropy.treesim.birth_death(birth_rate=1, death_rate=0.5, ntax=10)
    sim_traits(mle, verbose=True)
    mle.print_plot(display_width=70)

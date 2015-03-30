# import
import numpy as np
import dendropy 


def read_tree(treeFile):
    """Parsing a tree file that's in various possible formats
    Args:
    treeFile -- path to tree file
    Return:
    dendropy object
    """
    psblFormats = ['newick', 'nexus']
    for f in psblFormats:
        try:
            tree = dendropy.get_from_stream(treeFile, schema=f)
        except DataParseError:
            tree = None
            continue
        else:
            break
    
    if tree is None:
        msg = 'Do not recongnize tree file format (possible formats:{})'
        fmts = ','.join(psblFormats)
        raise DataParseError, msg.format(fmts)
    else:
        return tree



def sim_traits(tree, start=0, sigma=0.1, weight=0.5, verbose=False):
    """Trait simulation as detailed in:
    author = {Münkemüller, Tamara and Lavergne, Sébastien and Bzeznik,
              Bruno and Dray, Stéphane and Jombart,
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
            
    

if __name__ == '__main__':
    mle = dendropy.treesim.birth_death(birth_rate=1, death_rate=0.5, ntax=10)
    sim_traits(mle, verbose=True)
    mle.print_plot(display_width=70)

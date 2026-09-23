"""Build trees from distance matrices and read/write them in Newick format."""

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from skbio import DistanceMatrix, TreeNode
from skbio.tree import bme, nj, nni

LINKAGE_METHODS = ('average', 'ward', 'complete', 'single', 'weighted')


def hc_tree(df, method='average'):
    """
    Hierarchical clustering tree. "average" is UPGMA, the usual choice for Mash distances.
    Branch lengths are half the merge heights, so tip-to-tip distances match the clustering distances.

    :param df: square distance matrix (pandas DataFrame)
    :param method: scipy linkage method
    """
    condensed = squareform(df.to_numpy(), checks=False)
    linkage_matrix = linkage(condensed, method=method)
    return TreeNode.from_linkage_matrix(linkage_matrix, list(df.index))


def nj_tree(df):
    """Neighbour joining tree (unrooted). Negative branch lengths are set to 0."""
    return nj(DistanceMatrix(df.to_numpy(), list(df.index)))


def me_tree(df):
    """Balanced minimum evolution tree refined with nearest neighbour interchanges. Much faster than NJ."""
    dm = DistanceMatrix(df.to_numpy(), list(df.index))
    tree = bme(dm)
    # An unrooted tree with 3 tips has a single topology: nothing to refine (and scikit-bio's nni fails on it)
    return nni(tree, dm) if len(df) > 3 else tree


def read_newick(path):
    """
    Read a Newick tree. Underscores in names are kept as is.
    Numeric internal node labels are read as support values, so they are written back unquoted.
    """
    tree = TreeNode.read(str(path), format='newick', convert_underscores=False)
    tree.assign_supports()
    return tree


def quote(name):
    """Always single-quote labels so any character is allowed. Single quotes are escaped by doubling them."""
    return "'{}'".format(str(name).replace("'", "''"))


def format_length(length):
    return ':{:.10g}'.format(length) if length is not None and not np.isnan(length) else ''


def to_newick(tree):
    """
    Serialize a tree in Newick format.
    Iterative (no recursion limit on huge trees) and all labels are quoted, which scikit-bio's writer does not do.
    Internal nodes with a "support" value (bootstrap) get it as an unquoted label.
    """
    parts = list()
    # The stack holds nodes still to visit and text tokens to output once their children are done
    stack = [tree]
    while stack:
        item = stack.pop()
        if isinstance(item, str):
            parts.append(item)
            continue
        if item.children and getattr(item, 'support', None) is not None:
            label = str(item.support)  # Unquoted so tree viewers read it as a support value
        elif item.name is not None:
            label = quote(item.name)
        else:
            label = ''
        if item is not tree:
            label += format_length(item.length)
        if item.children:
            parts.append('(')
            stack.append(')' + label)
            for i, child in enumerate(reversed(item.children)):
                if i:
                    stack.append(',')
                stack.append(child)
        else:
            parts.append(label)
    return ''.join(parts) + ';\n'


def write_newick(tree, path):
    with open(path, 'w') as f:
        f.write(to_newick(tree))

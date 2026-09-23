"""
Bootstrap support for trees built from Mash distances.

There is no alignment to resample, so each replicate re-sketches the genomes with a different hash seed
("mash sketch -S"), which picks a different random subset of k-mers (the approach used by mashtree).
The support of a clade is the percentage of replicate trees that contain it.

Clades are stored as integer bitmasks over the tip names (bit i set = tip i is in the clade),
which keeps memory and comparisons cheap on trees with thousands of tips.
"""

from collections import Counter

ROOTED_TREES = ('hc',)  # Hierarchical clustering trees are rooted, NJ and ME trees are not


def clade_masks(tree, index):
    """
    :param index: dict {tip name: bit position}
    :return: dict {internal node: bitmask of the tips below it}, root excluded
    """
    masks = dict()
    for node in tree.postorder(include_self=False):
        if node.is_tip():
            masks[node] = 1 << index[node.name]
        else:
            mask = 0
            for child in node.children:
                mask |= masks[child]
            masks[node] = mask
    return {node: mask for node, mask in masks.items() if not node.is_tip()}


def canonical_split(mask, n_tips):
    """An unrooted split and its complement are the same bipartition. Keep the side without tip 0."""
    return mask ^ ((1 << n_tips) - 1) if mask & 1 else mask


def node_splits(tree, index, rooted):
    """
    :return: dict {internal node: clade (rooted) or bipartition (unrooted)}, trivial splits excluded
    """
    n_tips = len(index)
    splits = dict()
    for node, mask in clade_masks(tree, index).items():
        if not rooted:
            mask = canonical_split(mask, n_tips)
        size = mask.bit_count()
        if 1 < size < n_tips - (0 if rooted else 1):
            splits[node] = mask
    return splits


class SupportCounter:
    """Count how often the clades of a reference tree are found in replicate trees."""

    def __init__(self, tree, rooted):
        self.tree = tree
        self.rooted = rooted
        self.index = {name: i for i, name in enumerate(sorted(t.name for t in tree.tips()))}
        self.splits = node_splits(tree, self.index, rooted)
        self.wanted = set(self.splits.values())
        self.counts = Counter()
        self.replicates = 0

    def add(self, replicate_tree):
        self.add_splits(set(node_splits(replicate_tree, self.index, self.rooted).values()))

    def add_splits(self, found):
        """Add the clades (bitmasks) of one replicate tree, e.g. computed in another process."""
        self.counts.update(found & self.wanted)
        self.replicates += 1

    def assign(self):
        """Store the support (% of replicates) in the "support" attribute of each internal node."""
        for node, mask in self.splits.items():
            node.support = round(100 * self.counts[mask] / self.replicates) if self.replicates else None
        return self.tree

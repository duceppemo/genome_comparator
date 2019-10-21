#!/usr/local/env python3


__author__ = 'duceppemo'
__version__ = '0.1'


"""
https://www.biostars.org/p/97409/
"""

from ete3 import Tree
from argparse import ArgumentParser


class TreeCollapser(object):
    def __init__(self, args):
        # Arguments
        self.input_tree = args.input
        self.output_tree = args.output
        self.dist = args.distance

        # Run
        self.run()

    def run(self):
        # Parse tree file (Newick format)
        t = Tree(self.input_tree)

        # Now collapse nodes with distance smaller than
        TreeCollapser.collapse(t, self.dist)

        # Collapsed nodes are labeled, so you locate them and prune them
        for node in t.search_nodes(collapsed=True):
            for child in node.get_children():
                child.detach()

        # Write new collapsed tree to file
        t.write(format=1, outfile=self.output_tree)

    @staticmethod
    def mean(array):
        return sum(array)/float(len(array))

    @staticmethod
    def cache_distances(tree):
        """
        precalculate distances of all nodes to the root
        """
        node2rootdist = {tree: 0}
        for node in tree.iter_descendants('preorder'):
            node2rootdist[node] = node.dist + node2rootdist[node.up]
        return node2rootdist

    @staticmethod
    def collapse(tree, min_dist):
        # cache the tip content of each node to reduce the number of times the tree is traversed
        node2tips = tree.get_cached_content()
        root_distance = TreeCollapser.cache_distances(tree)

        for node in tree.get_descendants('preorder'):
            if not node.is_leaf():
                avg_distance_to_tips = TreeCollapser.mean([root_distance[tip]-root_distance[node]
                                                           for tip in node2tips[node]])

                if avg_distance_to_tips < min_dist:
                    # do whatever, ete support node annotation, deletion, labeling, etc.

                    # rename
                    # node.name += ' COLLAPSED avg_d:%g {%s}' %(avg_distance_to_tips,
                    #                                           ','.join([tip.name for tip in node2tips[node]]))
                    node.name += '%s {%s}' % ([tip.name for tip in node2tips[node]][0],
                                              ','.join([tip.name for tip in node2tips[node]][1:]))
                    # label
                    node.add_features(collapsed=True)

                    # set drawing attribute so they look collapsed when displayed with tree.show()
                    node.img_style['draw_descendants'] = False


if __name__ == '__main__':
    parser = ArgumentParser(description='Filter a blast xml file to keep the hits with a minimum % similarity'
                                        'Much quicker to do than rerun the whole blast')
    parser.add_argument('-i', '--input', metavar='sample.tree',
                        required=True,
                        help='Newick input tree')
    parser.add_argument('-o', '--output', metavar='sample_collapsed.tree',
                        required=True,
                        help='Newick output tree')
    parser.add_argument('-d', '--distance', metavar='0.01', type=float,
                        required=True,
                        help='Distance threshold. Nodes with distance smaller than this value will be collapsed.')

    # Get the arguments into an object
    arguments = parser.parse_args()

    TreeCollapser(arguments)

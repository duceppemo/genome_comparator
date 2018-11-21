#!/usr/local/bin python3


__author__ = 'duceppemo'
__version__ = '0.1'


import numpy as np
import pandas as pd
from os.path import basename
from time import time


class Dendro(object):

    def __init__(self, args):
        self.args = args
        self.input = args.input
        self.output = args.output
        self.nj_flag = args.nj

        # data
        self.df = pd.DataFrame()

        # Run all
        self.run()

    def run(self):
        start_time = time()
        ext = basename(self.input).split('.')[-1]

        # Parse input file
        t0 = time()
        print("Parsing input file to pandas dataframe...", end="", flush=True)
        if ext == 'xlsx' or ext == 'xls':
            self.df = pd.read_excel(self.input, index_col=0, header=0)
        elif ext == 'csv':
            self.df = pd.read_csv(self.input, index_col=0, header=0)
        elif ext == 'tsv':
            self.df = pd.read_csv(self.input, index_col=0, header=0, sep='\t')
        else:
            raise Exception("Invalid input file type."
                            "Input file should be in Excel ('.xlsx' or '.xls')"
                            "or text format ('.csv' or '.tsv')")

        print(" %s" % self.elapsed_time(time() - t0))

        # Check sample names for illegal characters for the tree file
        # self.fix_names(self.df)

        # Make tree from dataframe
        self.do_it(self.df)

        # Total time
        print("\nDone in %s" % self.elapsed_time(time() - start_time))

    def elapsed_time(self, seconds):
        """
        Transform a time value into a string
        :param seconds: Elapsed time in seconds
        :return: Formated time string
        """

        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)

        periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
        time_string = ''.join('{}{}'.format(int(np.round(value)), name) for name, value in periods if value)

        return time_string

    def fix_names(self, df):
        """
        Check sample names for illegal characters for the tree file
        :param df: pandas dataframe
        :return: fixed dataframe
        """

        import re

        t0 = time()
        print("Fixing illegal characters, if presents...", end="", flush=True)
        pattern = re.compile(r'[^0-9a-zA-Z_-]')
        new_headers = df.columns.str.replace(pattern, '')
        df.columns = new_headers  # Change column names
        df.index = new_headers  # Change row names
        print(" %s" % self.elapsed_time(time() - t0))

    def do_it(self, df):
        """
        Make tree from distance matrix
        :param df:
        :return:
        """

        # order data frame (square matix)
        print("Sorting dataframe...", end="", flush=True)
        df = df.reindex(sorted(df.columns), axis=1)  # columns
        df = df.reindex(sorted(df.columns), axis=0)  # rows

        # Get labels for tree and figure
        labels = df.columns.tolist()
        labels[:] = ['\'{}\''.format(x) for x in labels]  # Add single quotes around labels to allow special characters

        # Using scipy -> Takes about 8s to run with ~9,000 Salmonella genomes (fasta)
        self.make_hc_dendrogram(df, labels)

        # # Save tree as figure
        # t0 = time()
        # print("Saving tree as pdf figure...", end="", flush=True)
        # self.save_dendrogram_to_picture(linkage_matrix, labels)
        # print(" %s" % self.elapsed_time(time() - t0))

        # Using scikit-bio  -> Takes about 8h to run with ~9,000 Salmonella genomes (fasta)
        if self.nj_flag:
            self.make_nj_tree(df, labels)

    def make_hc_dendrogram(self, df, labels):
        from scipy.cluster.hierarchy import linkage, to_tree
        from scipy.spatial.distance import squareform

        # Convert square matrix to condensed distance matrix
        dists = squareform(df)

        # hierarchical clustering
        t0 = time()
        print("Making clusters...", end="", flush=True)
        linkage_matrix = linkage(dists, "ward")  # optimal_ordering=Ture is slow on large datasets
        print(" %s" % self.elapsed_time(time() - t0))

        # Create tree in Newick format
        t0 = time()
        print('Making dendrogram...', end="", flush=True)
        tree = to_tree(linkage_matrix, False)
        nw = self.getNewick(tree, "", tree.dist, labels)
        print(" %s" % self.elapsed_time(time() - t0))

        # save tree to file
        name = basename(self.input).split('.')[0]  # input file name without extension
        out_tree_file = self.output + '/' + name + '_hc.nwk'
        with open(out_tree_file, 'w') as out:
            out.write(nw)

    def save_dendrogram_to_picture(self, linkage_matrix, labels):
        from scipy.cluster.hierarchy import dendrogram
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(len(labels) / 500, len(labels) / 100))  # width, height, in inches

        # Create dendrogram
        dendrogram(linkage_matrix, labels=labels, orientation='left', leaf_font_size=2)

        # Output file name
        name = basename(self.input).split('.')[0]  # input file name without extension
        out_figure_file = self.output + '/' + name + '.pdf'

        # Customize the figure
        plt.title("Dendrogram")
        plt.tight_layout(rect=[0, 0, 0.99, 1])  # leave space on the right for leaf labels [left, bottom, right, top]

        fig.savefig(out_figure_file)

    def make_nj_tree(self, df, labels):
        from skbio import DistanceMatrix
        from skbio.tree import nj

        t0 = time()
        print("Converting dataframe to DistanceMatrix object...", end="", flush=True)
        dm = DistanceMatrix(df, labels)
        print(" %s" % self.elapsed_time(time() - t0))

        t0 = time()
        print("Creating nj tree...", end="", flush=True)
        nw = nj(dm, result_constructor=str)
        print(" %s" % self.elapsed_time(time() - t0))
        name = basename(self.input).split('.')[0]  # input file name without extension
        out_tree_file = self.output + '/' + name + '_nj.nwk'
        with open(out_tree_file, 'w') as out:
            out.write(nw)

    def getNewick(self, node, newick, parentdist, leaf_names):
        """
        https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
        :param node:
        :param newick:
        :param parentdist:
        :param leaf_names:
        :return:
        """
        if node.is_leaf():
            return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"

            newick = self.getNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = self.getNewick(node.get_right(), ",%s" % newick, node.dist, leaf_names)
            newick = "(%s" % newick

            return newick


if __name__ == "__main__":

    from argparse import ArgumentParser
    # https://docs.python.org/dev/library/argparse.html#action

    parser = ArgumentParser(description='Create dendrogram from a square distance matrix')
    parser.add_argument('-i', '--input', metavar='myfile.xlsx',
                        required=True,
                        help='Excel file (.xlsx) with the matrix.\n'
                             'First row and column have the sample names'
                             'First row and column shall all have the same names')
    parser.add_argument('-o', '--output', metavar='Out folder',
                        required=True,
                        help='Folder to hold the result files')
    # default=False is implied by action='store_true'
    parser.add_argument('-n', '--nj', action='store_true',
                        help='Output neighbour joining (nj) tree. Default "False".'
                             'Warning: can take several hours to complete on big matrices'
                             'e.g. it takes about 8h to run on a 9,000 x 9,000 matrix')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Dendro(arguments)
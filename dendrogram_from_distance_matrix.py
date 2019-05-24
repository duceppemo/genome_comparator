#!/usr/local/bin python3

import numpy as np
import pandas as pd
from os.path import basename
from time import time
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
from skbio import DistanceMatrix
from skbio.tree import nj
from sklearn import manifold
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt


class Dendro(object):

    def __init__(self, args):
        self.args = args
        self.input = args.input
        self.output = args.output
        self.nj_flag = args.nj
        self.pca_flag = args.pca

        # data
        self.df = pd.DataFrame()
        self.labels = list()

        # Run all
        self.run()

    def run(self):
        start_time = time()

        # Parse input file
        t0 = time()
        print("Parsing input file to pandas dataframe...", end="", flush=True)
        self.parse_input(self.input)
        nrow, ncol = self.df.shape
        print(" ({} x {}) {}".format(nrow, ncol, self.elapsed_time(time() - t0)))

        # Check sample names for illegal characters for the tree file
        # self.fix_names(self.df)

        # Sort dataframe
        print("Sorting dataframe...", end="", flush=True)
        t0 = time()
        self.sort_dataframe()
        print(" %s" % self.elapsed_time(time() - t0))

        # Get labels
        print("Getting labels...", end="", flush=True)
        t0 = time()
        self.get_labels()
        print(" %s" % self.elapsed_time(time() - t0))

        # Make tree from dataframe
        # Using scipy -> Takes about 8s to run with ~9,000 Salmonella genomes (fasta)
        print("Making hierarchical clustering tree...", end="", flush=True)
        t0 = time()
        self.make_hc_dendrogram(self.df, self.labels)
        print(" %s" % self.elapsed_time(time() - t0))

        # # Save tree as figure
        # t0 = time()
        # print("Saving tree as pdf figure...", end="", flush=True)
        # self.save_dendrogram_to_picture(linkage_matrix, labels)
        # print(" %s" % self.elapsed_time(time() - t0))

        # Make PCA
        if self.pca_flag:
            print("Making PCA...", end="", flush=True)
            t0 = time()
            self.make_pca(self.df)
            print(" %s" % self.elapsed_time(time() - t0))

        # Make tree from dataframe
        # Using scikit-bio  -> Takes about 8h to run with ~9,000 Salmonella genomes (fasta)
        if self.nj_flag:
            self.make_nj_tree(self.df, self.labels)

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

    def parse_input(self, input_file):
        """
        Parse input file into padans dataframe
        :param input_file: input file (square matrix)
        :return:
        """
        file_ext = basename(input_file).split('.')[-1]

        if file_ext == 'xlsx' or file_ext == 'xls':
            self.df = pd.read_excel(input_file, index_col=0, header=0)
        elif file_ext == 'csv':
            self.df = pd.read_csv(input_file, index_col=0, header=0)
        elif file_ext == 'tsv':
            self.df = pd.read_csv(input_file, index_col=0, header=0, sep='\t')
        else:
            raise Exception("Invalid input file type."
                            "Input file should be in Excel ('.xlsx' or '.xls')"
                            "or text format ('.csv' or '.tsv')")

    def sort_dataframe(self):
        """
        Make tree from distance matrix
        :return:
        """

        # order data frame (square matix)
        self.df = self.df.reindex(sorted(self.df.columns), axis='index')  # rows or index
        self.df = self.df.reindex(sorted(self.df.columns), axis='columns')  # columns

    def get_labels(self):
        """
        Get the column names from the pandas dataframe
        :return:
        """

        # Get labels for tree and figure
        self.labels = self.df.columns.tolist()

        # Add single quotes around labels to allow special characters
        self.labels[:] = ['\'{}\''.format(x) for x in self.labels]

    def make_hc_dendrogram(self, df, labels):
        """
        Make hierarchical clustering tree (newick format)
        :param df: pandas dataframe (square matrix)
        :param labels: column names (same as row names, same order too)
        :return:
        """

        # Convert square matrix to condensed matrix
        dists = squareform(df)

        # hierarchical clustering
        linkage_matrix = linkage(dists, "ward")  # optimal_ordering=Ture is slow on large datasets

        # Create tree in Newick format
        tree = to_tree(linkage_matrix, False)
        nw = self.get_newick(tree, "", tree.dist, labels)

        # save tree to file
        name = basename(self.input).split('.')[0]  # input file name without extension
        out_tree_file = self.output + '/' + name + '_hc.nwk'
        with open(out_tree_file, 'w') as out:
            out.write(nw)

    def save_dendrogram_to_picture(self, linkage_matrix, labels):
        """

        :param linkage_matrix:
        :param labels:
        :return:
        """

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
        """
        Make Neighbour Joining tree
        :param df: pandas dataframe (square matrix)
        :param labels: column names (same as row names, same order too)
        :return:
        """

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

    def get_newick(self, node, newick, parentdist, leaf_names):
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

            newick = self.get_newick(node.get_left(), newick, node.dist, leaf_names)
            newick = self.get_newick(node.get_right(), ",%s" % newick, node.dist, leaf_names)
            newick = "(%s" % newick

            return newick

    def make_pca(self, df):
        """

        :param df: pandas dataframe (square matrix)
        :return:
        """

        # dists = squareform(df)
        # dm = DistanceMatrix(df, labels)
        mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
        results = mds.fit(df)

        fig, ax = plt.subplots()
        coords = results.embedding_
        plt.subplots_adjust(bottom=0.1)
        plt.scatter(coords[:, 0], coords[:, 1], marker='o')

        species = [x.split('_')[1] for x in self.labels]
        for label, x, y in zip(species, coords[:, 0], coords[:, 1]):
            plt.annotate(
                label,
                xy=(x, y), xytext=(-10, 10),
                textcoords='offset points', ha='right', va='bottom', fontsize=6,
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        # texts = [plt.text(x, y, label,
        #                   ha='center', va='center') for label, x, y in zip(species, coords[:, 0], coords[:, 1])]
        # from adjustText import adjust_text
        # adjust_text(texts)

        # Output file name
        name = basename(self.input).split('.')[0]  # input file name without extension
        out_figure_file = self.output + '/' + name + '_PCA.png'

        fig.savefig(out_figure_file)


if __name__ == "__main__":

    from argparse import ArgumentParser
    # https://docs.python.org/dev/library/argparse.html#action

    parser = ArgumentParser(description='Create dendrogram from a square distance matrix')
    parser.add_argument('-i', '--input', metavar='my_square_matrix.tsv',
                        required=True,
                        help='Square matrix in tab-separated (.tsv) format.'
                             'First row and column have the sample names'
                             'First row and first column shall all have the same names in same order')
    parser.add_argument('-o', '--output', metavar='/output/folder/',
                        required=True,
                        help='Folder to hold the result files')
    # default=False is implied by action='store_true'
    parser.add_argument('--nj', action='store_true',
                        help='Output neighbour joining (nj) tree'
                             'Warning: can take several hours to complete on big matrices'
                             'e.g. it takes about 8h to run on a 9,000 x 9,000 matrix'
                             'Default "False"')
    parser.add_argument('--pca', action='store_true',
                        help='Output a PCA'
                             'Caution, experimental'
                             'Default "False"')

    # Get the arguments into an object
    arguments = parser.parse_args()

    Dendro(arguments)

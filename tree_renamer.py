from argparse import ArgumentParser


class TreeRenamer(object):
    def __init__(self, args):
        self.input = args.input
        self.output = args.output
        self.rename_table = args.rename_table

        # Data
        self.rename_dict = dict()

        # Do stuff
        self.parse_rename_talbe()
        self.rename_tree_file()

    def parse_rename_talbe(self):
        with open(self.rename_table, 'r') as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue
                old_name, new_name = line.split('\t')
                self.rename_dict[old_name] = new_name

    def rename_tree_file(self):
        with open(self.output, 'w') as out_fh:
            with open(self.input) as in_fh:
                for line in in_fh:
                    line = line.rstrip()
                    if not line:
                        continue
                    for old_name, new_name in self.rename_dict.items():
                        new_name = "'" + new_name + "'"
                        line = line.replace(old_name, new_name)
                    out_fh.write(line + '\n')


if __name__ == "__main__":
    parser = ArgumentParser(description='Create dendrogram from a square distance matrix')
    parser.add_argument('-i', '--input', metavar='/path_to/input_tree.nwk',
                        required=True, type=str,
                        help='Input tree in newick format')
    parser.add_argument('-o', '--output', metavar='/path_to/renamed_tree.nwk',
                        required=True, type=str,
                        help='Renamed tree in newick format')
    parser.add_argument('-r', '--rename-table', metavar='/path_to/rename_table.tsv',
                        required=True, type=str,
                        help='Renaming talbe. File with two columns separated with "tab" character. '
                             'First column is actual sample names in tree and second is new names.')
    # Get the arguments into an object
    arguments = parser.parse_args()

    TreeRenamer(arguments)

#!/usr/local/bin python3

import os
import pathlib
import subprocess
from argparse import ArgumentParser
from concurrent import futures
from glob import glob
from multiprocessing import cpu_count


class MashPhylo(object):

    def __init__(self, args):

        self.args = args
        self.input = args.input
        self.output = args.output
        self.cpu = args.threads
        self.sketch_size = args.sketch_size
        self.kmer_size = args.kmer_size

        # Performance
        mcpu = cpu_count()
        if self.cpu > mcpu:
            self.cpu = mcpu

        # Data
        self.input_file_list = list()
        self.sketch_list = list()
        self.distance_list = list()

        # Look for input sequence files recursively
        accepted_extentions = ['.fna', '.fna.gz',
                               '.fasta', '.fasta.gz',
                               '.fa', '.fa.gz',
                               '.fastq', '.fastq.gz',
                               '.fq', 'fq.gz']

        for root, directories, filenames in os.walk(self.input):
            for filename in filenames:
                if filename.endswith(tuple(accepted_extentions)):  # accept a tuple or string
                    self.input_file_list.append(os.path.join(root, filename))

        # Run
        self.run()

    def run(self):
        self.checks()
        self.create_folders(self.output)
        self.parallel_sketch_it(self.input_file_list)
        self.paste_sketches(os.path.join(self.output, 'sketches'))
        self.parallel_dist_it(self.sketch_list)
        self.create_distance_matrix(self.distance_list)
        self.create_tree(os.path.join(self.output, 'all_dist.tsv'))

    def checks(self):
        if not os.path.isdir(self.input):
            raise Exception('An input folder is required')
        if not os.path.exists(self.input):
            raise Exception('The provided input folder does not exists')

        # Check self.sketch_size

        # Check self.kmer_size

        # Check if input seqence files were found
        if not self.input_file_list:
            raise Exception('No valid sequence files detected in provided input folder')

        # Need a least 3 samples to make a tree
        if len(self.input_file_list) < 3:
            raise Exception('At least 3 samples are required to build a tree')

    def create_folders(self, output):
        # Create output folder is it does not exist
        pathlib.Path(output).mkdir(parents=True, exist_ok=True)
        # Create output folder for mash sketches
        pathlib.Path(os.path.join(output, 'sketches')).mkdir(parents=True, exist_ok=True)
        # Create output folder for mash distances
        pathlib.Path(os.path.join(output, 'distances')).mkdir(parents=True, exist_ok=True)
        # Create output folder for dendrogram
        pathlib.Path(os.path.join(output, 'tree')).mkdir(parents=True, exist_ok=True)

    def sketch_it(self, sequence_file):
        sample_name = os.path.basename(sequence_file).split('.')[0]
        output_file = os.path.join(self.output, 'sketches', sample_name)

        # Use one CPU per process
        fasta_ext = ['.fna', '.fna.gz', '.fasta', '.fasta.gz', '.fa', '.fa.gz']
        fastq_ext = ['.fastq', '.fastq.gz', '.fq', 'fq.gz']

        if sequence_file.endswith(tuple(fasta_ext)):
            subprocess.run(['mash', 'sketch',
                            '-k', str(self.kmer_size),
                            '-s', str(self.sketch_size),
                            '-p', '1',
                            '-o', output_file,
                            sequence_file])
        elif sequence_file.endswith(tuple(fastq_ext)):
            subprocess.run(['mash', 'sketch',
                            '-k', str(self.kmer_size),
                            '-s', str(self.sketch_size),
                            '-p', '1',
                            '-r',
                            '-m', '2',
                            '-o', output_file,
                            sequence_file])
        else:
            pass  # Have been taking care of earlier

    def parallel_sketch_it(self, file_list):
        """

        :param file_list: list of files (absolute paths)
        :return:
        """
        # we would want to use the ProcessPoolExecutor for CPU intensive tasks.
        # The ThreadPoolExecutor is better suited for network operations or I/O.
        with futures.ThreadPoolExecutor(max_workers=self.cpu) as executor:
            future_to_sketch = {executor.submit(self.sketch_it(s)): s for s in file_list}
            # future_to_sketch = {executor.submit(self.sketch_it(s)): s for s in file_list}
            # for future in futures.as_completed(future_to_sketch):
            #     sketch = future_to_sketch[future]

    def paste_sketches(self, sketch_folder):
        search_pattern = sketch_folder + '/*.msh'
        self.sketch_list = glob(search_pattern)
        # write list to file
        with open(os.path.join(self.output, 'sketch.list'), 'w') as f:
            for item in self.sketch_list:
                f.write('{}\n'.format(item))

        proc = subprocess.run(['mash', 'paste',
                               os.path.join(self.output, 'all.msh'),
                               '-l', os.path.join(self.output, 'sketch.list')])

    def dist_it(self, sketch_file):
        sample_name = os.path.basename(sketch_file).split('.')[0]

        proc = subprocess.Popen(['mash', 'dist',
                                 '-p', '1',
                                 '-t',
                                 os.path.join(self.output, 'all.msh'),
                                 sketch_file], stdout=subprocess.PIPE)

        stdout, stderr = proc.communicate()

        # remove path and file extension from file names
        my_sample = stdout.split(b'\t')[1]
        folder = os.path.dirname(my_sample) + b'/'
        ext = b'.' + b'.'.join(my_sample.split(b'.')[1:])  # Account for double extension like ".fastq.gz"
        stdout = stdout.replace(folder, b'').replace(ext, b'')

        # TODO -> parse into PANDA dataframe instead of writting to file. Maybe?
        # Write distance file
        output_file = os.path.join(self.output, 'distances', sample_name) + '_dist.tsv'
        self.distance_list.append(output_file)
        with open(output_file, 'w') as f:
            f.write(stdout.decode('ascii'))

    def parallel_dist_it(self, sketch_list):
        with futures.ThreadPoolExecutor(max_workers=self.cpu) as executor:
            future_to_dist = {executor.submit(self.dist_it(d)): d for d in sketch_list}
            # for future in futures.as_completed(future_to_dist):
            #     sketch = future_to_dist[future]

    def create_distance_matrix(self, distance_list):
        header = open(distance_list[0], 'r').readline()

        with open(os.path.join(self.output, 'all_dist.tsv'), 'w') as outfile:
            # write header
            outfile.write(header)
            for d in distance_list:
                with open(d, 'r') as infile:
                    for line in infile:
                        if not line.startswith('#'):
                            outfile.write(line)

    def create_tree(self, matrix_file):
        subprocess.run(['python3', 'dendrogram_from_distance_matrix.py',
                        '-i', matrix_file,
                        '-o', os.path.join(self.output, 'tree'),
                        '--nj'])


if __name__ == "__main__":
    # https://docs.python.org/dev/library/argparse.html#action
    max_cpu = cpu_count()

    parser = ArgumentParser(description='Create dendrogram from a square distance matrix')
    parser.add_argument('-i', '--input', metavar='/input/folder',
                        required=True,
                        help='Input folder containing the fasta or fastq files')
    parser.add_argument('-o', '--output', metavar='/output/folder',
                        required=True,
                        help='Folder to hold the result files')
    parser.add_argument('-t', '--threads', metavar=max_cpu,
                        required=False,
                        type=int, default=max_cpu,
                        help='Number of CPU.'
                             'Default is maximum CPU available({})'.format(max_cpu))
    parser.add_argument('-k', '--kmer_size',
                        required=False,
                        type=int, default=21,
                        help='kmer size to use for Mash.'
                             'Default 21.')
    parser.add_argument('-s', '--sketch_size',
                        required=False,
                        type=int, default=10000,
                        help='Sketch size to use for Mash'
                             'Default 10,000.')
    # default=False is implied by action='store_true'
    parser.add_argument('--nj', action='store_true',
                        help='Output neighbour joining (nj) tree. Default "False".'
                             'Warning: can take several hours to complete on big matrices'
                             'e.g. it takes about 8h to run on a 9,000 x 9,000 matrix')
    parser.add_argument('--pca', action='store_true',
                        help='Output PCA. Default "False".'
                             'Not very useful when too many samples.'
                             'Still experimental')

    # Get the arguments into an object
    arguments = parser.parse_args()

    MashPhylo(arguments)

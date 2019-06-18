#!/usr/local/bin python3

import os
import pathlib
import subprocess
from argparse import ArgumentParser
from concurrent import futures
from multiprocessing import cpu_count
from time import time
import numpy as np
import gzip


# TODO -> use logger
# TODO -> Maybe add progress bar for the sketching? # https://renzo.lucioni.xyz/parallel-progress/


class SampleObject(object):
    def __init__(self, sample_name, file_type, file_path, sketch_file, dist_file, est_size, est_cov, contigs):
        # Create seq object with its attributes
        self.sample_name = sample_name
        self.file_type = file_type  # fastq or fasta
        self.file_path = [file_path]  # list
        self.sketch_file = sketch_file
        self.dist_file = dist_file
        self.est_size = est_size
        self.est_cov = est_cov
        self.contigs = contigs


class MashPhylo(object):

    def __init__(self, args):

        # Paths
        self.home = str(pathlib.Path.home())
        self.cwd = os.getcwd()
        self.script_path = os.path.dirname(os.path.abspath(__file__))

        # Command line's arguments
        self.args = args
        self.input = self.fix_paths(args.input)
        self.output = self.fix_paths(args.output)
        self.cpu = args.threads
        self.sketch_size = args.sketch_size
        self.kmer_size = args.kmer_size
        self.nj = args.nj
        self.pca = args.pca

        # Performance
        mcpu = cpu_count()
        if self.cpu > mcpu:
            self.cpu = mcpu

        # Data
        self.input_dict = dict()

        # Run
        self.run()

    def run(self):
        # To log total runtime
        start_time = time()

        my_dict = self.input_dict
        self.checks()
        self.get_samples(self.input)
        MashPhylo.create_output_folders(self.output)

        print("Sketching files...", end="", flush=True)
        t0 = time()
        self.parallel_sketch_it(my_dict)
        print(" %s" % self.elapsed_time(time() - t0))

        MashPhylo.write_stats(self.input_dict, os.path.join(self.output, 'sample_stats.txt'))
        self.paste_sketches(my_dict)

        print("Pasting all sketches together...", end="", flush=True)
        t0 = time()
        self.parallel_dist_it(my_dict)
        print(" %s" % self.elapsed_time(time() - t0))

        print("Computing distances...", end="", flush=True)
        t0 = time()
        self.create_distance_matrix(my_dict)
        print(" %s" % self.elapsed_time(time() - t0))

        self.clean_samples_names(os.path.join(self.output, 'all_dist.tsv'))

        print("Making tree(s)...")
        self.create_tree(os.path.join(self.output, 'all_dist.tsv'))

        # Total time
        print("\nDone in %s" % MashPhylo.elapsed_time(time() - start_time))

    def checks(self):
        if not os.path.isdir(self.input):
            raise Exception('The provided input path is not a directory')
        if not os.path.exists(self.input):
            raise Exception('The provided input directory does not exists')

        # Check self.sketch_size
        # Make sure that there are not too many sketches.
        # For example it doesn't make sense to split a 10,000bp genome into 10,000 sketches!
        # Also you can't create more sketches that the total bp!
        # Maybe compare it to file size?

        # Check self.kmer_size is not out of allowed sized by Mash

    def fix_paths(self, in_path):
        """
        Convert relative paths to absolute path. Also convert tild ("~") to its absolute value.

        :param in_path: string; path
        :return: out_path; string; absolute path of in_path
        """
        out_path = in_path

        if os.path.dirname(in_path).startswith('./'):
            out_path = in_path.replace('.', self.cwd, 1)  # only replace the first one (not in the file extension)
        elif os.path.dirname(in_path).startswith('../'):
            out_path = os.path.dirname(in_path.replace('..', self.cwd, 1))
        elif in_path.startswith('~'):
            out_path = in_path.replace('~', self.home)

        return out_path

    def get_samples(self, input_folder):
        """
        Creates a dictionary of 'sample objects'.

        :param input_folder: string; absolute path to input folder containing fastq and/or fasta files
        :return:
        """
        accepted_extentions = ['.fna', '.fna.gz',
                               '.fasta', '.fasta.gz',
                               '.fa', '.fa.gz',
                               '.fastq', '.fastq.gz',
                               '.fq', '.fq.gz']

        # Look for input sequence files recursively
        for root, directories, filenames in os.walk(input_folder):
            for filename in filenames:
                if filename.endswith(tuple(accepted_extentions)):  # accept a tuple or string
                    sample_name = os.path.basename(filename).split('.')[0].split('_')[0]
                    file_path = os.path.join(root, filename)
                    file_type = os.path.basename(filename).split('.')[1]
                    if '.' + file_type not in accepted_extentions:
                        raise Exception('Don\'t use dot "." in your file names')

                    sample_object = SampleObject(sample_name, file_type, file_path, None, None, 'N/A', 'N/A', None)

                    if sample_name not in self.input_dict.keys():
                        self.input_dict[sample_name] = sample_object
                    # multiple files per samples are allowed (e.g. R1 and R2)
                    elif sample_name in self.input_dict.keys():
                        if file_type == self.input_dict[sample_name].file_type:
                            self.input_dict[sample_name].file_path.append(file_path)
                        else:
                            raise Exception('Sample {} has both fastq and fasta files. Keep only one file type'.format(
                                sample_name))

        # Check if input sequence files were found
        if not self.input_dict:
            raise Exception('No valid sequence files detected in provided input folder')

        # Need a least 3 samples to make a tree
        if len(self.input_dict.keys()) < 3:
            raise Exception('At least 3 samples are required to build a tree')

    @staticmethod
    def elapsed_time(seconds):
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

    @staticmethod
    def create_output_folders(output):
        """
        Create the output folder and subfolders

        :param output: string; absolute path to output folder
        :return:
        """
        # Create output folder is it does not exist
        pathlib.Path(output).mkdir(parents=True, exist_ok=True)
        # Create output folder for mash sketches
        pathlib.Path(os.path.join(output, 'sketches')).mkdir(parents=True, exist_ok=True)
        # Create output folder for mash distances
        pathlib.Path(os.path.join(output, 'distances')).mkdir(parents=True, exist_ok=True)
        # Create output folder for dendrogram
        pathlib.Path(os.path.join(output, 'tree')).mkdir(parents=True, exist_ok=True)


    @staticmethod
    def get_assembly_info(fasta_file):
        """
        Compute assembly size in bp and number of contigs
        :param fasta_file: string; absolute path to assembly file
        :return: string; string representation of the assembly size in bp
                 string; string representation of number of contigs in the assembly
        """

        assembly_size = 0
        contigs = 0
        with gzip.open(fasta_file, 'r') if fasta_file.endswith('.gz') else open(fasta_file, 'r') as f:
            for line in f:
                line = line.rstrip()
                if len(line) == 0:  # skip empty lines
                    continue

                if line.startswith('>'):
                    contigs += 1
                else:
                    assembly_size += len(line)
        return str(assembly_size), str(contigs)

    def sketch_it(self, info):
        output_file = os.path.join(self.output, 'sketches', info.sample_name)
        info.sketch_file = output_file + '.msh'

        fasta_ext = ['fna', 'fasta', 'fa']
        fastq_ext = ['fastq', 'fq']
        is_fastq = False
        cmd = list()

        if info.file_type in fasta_ext:
            # Use one CPU per process (-p 1)
            cmd = ['mash', 'sketch',
                   '-k', str(self.kmer_size),
                   '-s', str(self.sketch_size),
                   '-p', '1',
                   '-o', output_file] + info.file_path
        elif info.file_type in fastq_ext:
            is_fastq = True
            # The option p will be ignored with r
            cmd = ['mash', 'sketch',
                   '-k', str(self.kmer_size),
                   '-s', str(self.sketch_size),
                   '-r',
                   '-m', '2',
                   '-I', info.sample_name,
                   '-o', output_file] + info.file_path
        else:
            pass  # Have been taking care of earlier

        # subprocess.run(cmd)  # Makes too much info on screen
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout, stderr = p1.communicate()

        # if fastq, output "Estimated genome size" and "Estimated coverage" into log file
        if is_fastq:
            stats = stdout.split(b'\n')
            info.est_size = str(int(float(stats[0].split(b' ')[-1].decode('ascii'))))
            info.est_cov = stats[1].split(b' ')[-1].decode('ascii')
        else:
            info.est_size, info.contigs = MashPhylo.get_assembly_info(info.file_path[0])

    def parallel_sketch_it(self, sample_dict):
        # we would want to use the ProcessPoolExecutor for CPU intensive tasks.
        # The ThreadPoolExecutor is better suited for network operations or I/O.

        with futures.ThreadPoolExecutor(max_workers=self.cpu) as executor:
            args = (info for sample, info in sample_dict.items())
            for results in executor.map(self.sketch_it, args):
                pass

    @staticmethod
    def write_stats(sample_dict, stats_file):
        with open(stats_file, 'w') as f:
            for name, info in sample_dict.items():
                if info.file_type == 'fastq':
                    pair = ''
                    if len(info.file_path) == 2:
                        pair = 'PE'
                    else:
                        pair = 'SE'
                    f.write("{} ({} {}):\n\tEstimated genome size: {}\n\tEstimated coverage: {}\n\n".format(
                        info.sample_name, pair, info.file_type, info.est_size, info.est_cov))
                else:  # fasta
                    f.write("{} ({}):\n\tAssembly size: {}\n\tNumber of contigs: {}\n\n".format(
                        info.sample_name, info.file_type, info.est_size, info.contigs))

    def paste_sketches(self, sample_dict):
        list_file = os.path.join(self.output, 'sketch.list')
        all_msh = os.path.join(self.output, 'all.msh')

        # write list to file
        with open(list_file, 'w') as f:
            for sample, info in sample_dict.items():
                f.write('{}\n'.format(info.sketch_file))

        proc = subprocess.run(['mash', 'paste',
                               all_msh,
                               '-l', list_file])

    def dist_it(self, info):
        info.dist_file = os.path.join(self.output, 'distances', info.sample_name) + '_dist.tsv'

        proc = subprocess.Popen(['mash', 'dist',
                                 '-p', '1',
                                 '-t',
                                 os.path.join(self.output, 'all.msh'),
                                 info.sketch_file], stdout=subprocess.PIPE)

        # run the process and caption standard output and standard error
        stdout, stderr = proc.communicate()

        # TODO -> parse into PANDA data frame instead of writing to file. Maybe?
        # Write distance file
        with open(info.dist_file, 'w') as f:
            f.write(stdout.decode('ascii'))

    def parallel_dist_it(self, sample_dict):
        with futures.ThreadPoolExecutor(max_workers=self.cpu) as executor:
            args = (info for sample, info in sample_dict.items())
            for results in executor.map(self.dist_it, args):
                pass

    def create_distance_matrix(self, sample_dict):
        distance_list = list()
        for sample, info in sample_dict.items():
            distance_list.append(info.dist_file)
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
        cmd = ['python3', self.script_path + '/dendrogram_from_distance_matrix.py',
               '-i', matrix_file,
               '-o', os.path.join(self.output, 'tree')]
        if self.nj:
            cmd.append('--nj')
        if self.pca:
            cmd.append('--pca')

        subprocess.run(cmd)

    def clean_samples_names(self, infile):
        # Create a new dictionary for renaming
        rename_dict = dict()
        for sample, info in self.input_dict.items():
            for pa in info.file_path:
                rename_dict[pa] = sample

        with open(infile + '.tmp', 'w') as tmpfile:
            with open(infile, 'r') as f:
                for line in f:
                    for p, n in rename_dict.items():
                        line = line.replace(p, n)
                    tmpfile.write(line)

        # rename file
        os.rename(infile + '.tmp', infile)


if __name__ == "__main__":
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

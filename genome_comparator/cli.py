"""Command line entry points."""

import logging
import os
import sys
from argparse import ArgumentParser, ArgumentTypeError, ArgumentDefaultsHelpFormatter
from pathlib import Path

from . import __version__, matrix, ordination, trees
from .mash import MAX_KMER_SIZE, MashError
from .matrix import MatrixError
from .pipeline import MIN_SAMPLES, MashPhylo, analyze_matrix, step
from .samples import SampleError
from .tree_tools import collapse, read_rename_table, rename_tips

log = logging.getLogger('genome_comparator')

# Expected errors are reported without a traceback
USER_ERRORS = (MashError, MatrixError, SampleError, ValueError, OSError)


def setup_logging(verbose=False, log_file=None):
    """Configure the package logger only, so importing genome_comparator never alters the root logger."""
    for handler in list(log.handlers):
        log.removeHandler(handler)
        handler.close()
    handlers = [logging.StreamHandler(sys.stderr)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode='w'))
    formatter = logging.Formatter('%(asctime)s %(levelname)-7s %(message)s', datefmt='%H:%M:%S')
    for handler in handlers:
        handler.setFormatter(formatter)
        log.addHandler(handler)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)


def int_range(low, high=None):
    def check(value):
        try:
            value = int(value)
        except ValueError:
            raise ArgumentTypeError('"{}" is not an integer'.format(value))
        if value < low or (high is not None and value > high):
            raise ArgumentTypeError('must be between {} and {}'.format(low, high) if high is not None
                                    else 'must be >= {}'.format(low))
        return value
    return check


def available_cpus():
    try:
        return len(os.sched_getaffinity(0))  # Respects cgroup/taskset limits on Linux (e.g. SLURM jobs)
    except AttributeError:
        return os.cpu_count() or 1


def add_tree_arguments(parser):
    group = parser.add_argument_group('trees and ordination')
    group.add_argument('--linkage', choices=trees.LINKAGE_METHODS, default='average',
                       help='Hierarchical clustering method for the "_hc" tree. "average" is UPGMA.')
    group.add_argument('--nj', action='store_true',
                       help='Also build a neighbour joining tree. Slower than --me on very large datasets.')
    group.add_argument('--me', action='store_true',
                       help='Also build a balanced minimum evolution tree (with NNI). '
                            'Similar to NJ but much faster on large datasets.')
    group.add_argument('--pcoa', '--pca', dest='pcoa', action='store_true',
                       help='Also run a principal coordinates analysis (PCoA) and save an interactive html plot.')
    group.add_argument('--metadata', metavar='metadata.tsv',
                       help='Tab-separated file, first column is the sample name. '
                            'Extra columns are shown when hovering over PCoA points.')
    group.add_argument('--color-by', metavar='COLUMN',
                       help='Metadata column used to colour the PCoA points.')


def check_tree_arguments(parser, args):
    if args.color_by and not args.metadata:
        parser.error('--color-by requires --metadata')
    if args.metadata and not args.pcoa:
        parser.error('--metadata is only used with --pcoa')
    if args.metadata:
        try:
            ordination.read_metadata(args.metadata, args.color_by)  # Fail now rather than after hours of work
        except (ValueError, OSError) as e:
            parser.error(str(e))


def tree_kwargs(args):
    return dict(linkage=args.linkage, nj=args.nj, me=args.me, pcoa=args.pcoa,
                metadata=args.metadata, color_by=args.color_by)


def run_safely(func):
    try:
        func()
    except USER_ERRORS as e:
        log.error('%s', e)
        log.debug('Traceback:', exc_info=True)  # Shown with --verbose
        sys.exit(1)
    except KeyboardInterrupt:
        log.error('Interrupted')
        sys.exit(130)


def main(argv=None):
    """genome-comparator: compare genomes from a folder of fasta/fastq files."""
    max_cpu = available_cpus()
    parser = ArgumentParser(prog='genome-comparator', formatter_class=ArgumentDefaultsHelpFormatter,
                            description='Compare genomes (assemblies or reads) with Mash and build '
                                        'a distance matrix, trees and an optional PCoA plot.')
    parser.add_argument('-i', '--input', metavar='/input/folder', required=True,
                        help='Folder containing the fasta or fastq files (searched recursively)')
    parser.add_argument('-o', '--output', metavar='/output/folder', required=True,
                        help='Folder to hold the result files')
    parser.add_argument('-t', '--threads', metavar='N', type=int_range(1), default=max_cpu,
                        help='Number of threads')
    parser.add_argument('-k', '--kmer-size', '--kmer_size', dest='kmer_size', type=int_range(1, MAX_KMER_SIZE),
                        default=21, help='k-mer size used by Mash')
    parser.add_argument('-s', '--sketch-size', '--sketch_size', dest='sketch_size', type=int_range(1),
                        default=10000, help='Number of min-hashes per sketch')
    parser.add_argument('-m', '--min-copies', type=int_range(1), default=2,
                        help='Reads only: minimum copies of a k-mer to be included in the sketch '
                             '(filters out sequencing errors)')
    add_tree_arguments(parser)
    parser.add_argument('-b', '--bootstrap', metavar='N', type=int_range(0), default=0,
                        help='Number of bootstrap replicates for tree support values. Each replicate sketches all '
                             'the samples again with a different hash seed, so N replicates take about N times '
                             'longer than a normal run.')
    parser.add_argument('--phylip', action='store_true',
                        help='Also save the distance matrix in phylip format (for rapidnj, fastme, etc.)')
    parser.add_argument('--force', action='store_true',
                        help='Sketch all samples again, even if up-to-date sketches exist in the output folder')
    parser.add_argument('--clean', action='store_true',
                        help='Remove the individual sketch files at the end')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show debug messages')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args(argv)
    check_tree_arguments(parser, args)

    if args.threads > max_cpu:
        args.threads = max_cpu

    output = Path(args.output).expanduser()
    output.mkdir(parents=True, exist_ok=True)
    setup_logging(args.verbose, output / 'genome_comparator.log')
    log.info('genome_comparator %s: %s', __version__, ' '.join(sys.argv))

    run_safely(MashPhylo(args.input, args.output, threads=args.threads, kmer_size=args.kmer_size,
                         sketch_size=args.sketch_size, min_copies=args.min_copies, phylip=args.phylip,
                         force=args.force, clean=args.clean, bootstrap=args.bootstrap,
                         **tree_kwargs(args)).run)


def mash_phylo_main(argv=None):
    """Deprecated name of the main command, kept so existing scripts keep working."""
    sys.stderr.write('Warning: "mash-phylo" is deprecated and will be removed in a future version. '
                     'Use "genome-comparator" instead (same options).\n')
    main(argv)


def dendrogram_main(argv=None):
    """dendrogram-from-matrix: build trees / PCoA from an existing square distance matrix."""
    parser = ArgumentParser(prog='dendrogram-from-matrix', formatter_class=ArgumentDefaultsHelpFormatter,
                            description='Build trees and an optional PCoA plot from a square distance matrix. '
                                        'First row and first column hold the sample names.')
    parser.add_argument('-i', '--input', metavar='my_square_matrix.tsv', required=True,
                        help='Square distance matrix (.tsv, .csv, .xlsx or .xls)')
    parser.add_argument('-o', '--output', metavar='/output/folder', required=True,
                        help='Folder to hold the result files')
    add_tree_arguments(parser)
    parser.add_argument('-v', '--verbose', action='store_true', help='Show debug messages')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args(argv)
    check_tree_arguments(parser, args)
    setup_logging(args.verbose)

    def run():
        with step('Reading distance matrix'):
            df = matrix.read_matrix(args.input)
        log.info('%d x %d matrix', *df.shape)
        if len(df) < MIN_SAMPLES:
            raise MatrixError('At least {} samples are required to build a tree'.format(MIN_SAMPLES))
        name = Path(args.input).stem
        for out in analyze_matrix(df, Path(args.output).expanduser(), name, **tree_kwargs(args)):
            log.info('Output: %s', out)

    run_safely(run)


def collapse_main(argv=None):
    """tree-collapser: collapse clades whose tips are closer than a distance threshold."""
    parser = ArgumentParser(prog='tree-collapser',
                            description='Collapse clades whose average distance to their tips is smaller '
                                        'than a threshold. Collapsed clades are replaced by a single tip '
                                        'named "<first tip> {<other tips>}".')
    parser.add_argument('-i', '--input', metavar='tree.nwk', required=True, help='Newick input tree')
    parser.add_argument('-o', '--output', metavar='tree_collapsed.nwk', required=True, help='Newick output tree')
    parser.add_argument('-d', '--distance', metavar='0.01', type=float, required=True,
                        help='Distance threshold. Clades with an average distance to their tips smaller than '
                             'this value are collapsed.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args(argv)
    setup_logging()

    def run():
        tree = trees.read_newick(args.input)
        n = collapse(tree, args.distance)
        trees.write_newick(tree, args.output)
        log.info('Collapsed %d clade(s)', n)

    run_safely(run)


def rename_main(argv=None):
    """tree-renamer: rename tree tips from a two-column table."""
    parser = ArgumentParser(prog='tree-renamer', description='Rename the tips of a Newick tree.')
    parser.add_argument('-i', '--input', metavar='input_tree.nwk', required=True, help='Input tree in Newick format')
    parser.add_argument('-o', '--output', metavar='renamed_tree.nwk', required=True,
                        help='Renamed tree in Newick format')
    parser.add_argument('-r', '--rename-table', metavar='rename_table.tsv', required=True,
                        help='Tab-separated file with two columns: current tip name, new name')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args(argv)
    setup_logging()

    def run():
        tree = trees.read_newick(args.input)
        not_found = rename_tips(tree, read_rename_table(args.rename_table))
        trees.write_newick(tree, args.output)
        if not_found:
            log.warning('%d name(s) from the rename table were not found in the tree: %s',
                        len(not_found), ', '.join(sorted(not_found)[:10]))

    run_safely(run)


if __name__ == '__main__':
    main()

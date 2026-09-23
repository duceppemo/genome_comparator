"""The genome comparison pipeline: sketch -> paste -> pairwise distances -> trees / PCoA."""

import csv
import logging
import multiprocessing
import tempfile
from collections import deque
from concurrent import futures
from contextlib import contextmanager
from pathlib import Path
from time import time

import pandas as pd

from . import mash, matrix, ordination, trees
from .bootstrap import ROOTED_TREES, SupportCounter, node_splits
from .samples import SampleError, find_samples

log = logging.getLogger(__name__)

MIN_SAMPLES = 3  # Minimum number of samples to build a tree
MAX_TREE_WORKERS = 8  # Processes building bootstrap trees in parallel


def elapsed_time(seconds):
    """Format a duration, e.g. "1h2m3s"."""
    minutes, seconds = divmod(round(seconds), 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    periods = [('d', days), ('h', hours), ('m', minutes), ('s', seconds)]
    return ''.join('{}{}'.format(value, name) for name, value in periods if value) or '0s'


@contextmanager
def step(message):
    """Log the start of a step and how long it took."""
    log.info('%s...', message)
    t0 = time()
    yield
    log.info('%s done in %s', message, elapsed_time(time() - t0))


@contextmanager
def thread_pool(threads):
    """
    ThreadPoolExecutor that drops its queued jobs on Ctrl-C. By default, the executor still runs every job
    submitted before the interruption, so stopping a run on thousands of samples could take as long as the run.
    """
    with futures.ThreadPoolExecutor(max_workers=threads) as executor:
        try:
            yield executor
        except KeyboardInterrupt:
            executor.shutdown(cancel_futures=True)
            raise


TREE_BUILDERS = {
    'hc': ('hierarchical clustering tree', lambda df, linkage: trees.hc_tree(df, linkage)),
    'me': ('minimum evolution tree', lambda df, linkage: trees.me_tree(df)),
    'nj': ('neighbour joining tree', lambda df, linkage: trees.nj_tree(df)),
}


def build_trees(df, linkage='average', nj=False, me=False, verbose=True):
    """:return: dict {tree kind: TreeNode}. The "hc" tree is always built."""
    kinds = ['hc'] + (['me'] if me else []) + (['nj'] if nj else [])
    built = dict()
    for kind in kinds:
        label, builder = TREE_BUILDERS[kind]
        if verbose:
            with step('Building {}{}'.format(label, ' ({})'.format(linkage) if kind == 'hc' else '')):
                built[kind] = builder(df, linkage)
        else:
            built[kind] = builder(df, linkage)
    return built


def worker_context():
    """
    Never "fork" worker processes: forking a process that runs threads (thread pools, the process pool's own
    management thread, test runners...) can deadlock the child. "fork" is the Linux default before Python 3.14.
    """
    method = 'forkserver' if 'forkserver' in multiprocessing.get_all_start_methods() else 'spawn'
    return multiprocessing.get_context(method)


def replicate_splits(rep_df, index, rooted, linkage, nj, me):
    """
    Build the trees of one bootstrap replicate and return their clades as bitmasks.
    Runs in a worker process: only the small sets of bitmasks are sent back, not the trees.

    :param rooted: dict {tree kind: True if the tree is rooted}
    """
    built = build_trees(rep_df, linkage, nj, me, verbose=False)
    return {kind: set(node_splits(tree, index, rooted[kind]).values()) for kind, tree in built.items()}


def add_bootstrap_support(built, replicates, linkage='average', nj=False, me=False, workers=1):
    """
    Build the same trees from each replicate distance matrix and set the support of every clade of the
    reference trees to the % of replicate trees containing it.

    Replicate trees are built in parallel worker processes while the next replicate matrices are produced.
    At most workers + 1 replicate matrices are kept in memory at once.

    :param built: reference trees from build_trees()
    :param replicates: iterable of replicate distance matrices (same samples as the reference)
    :param workers: number of processes building replicate trees
    """
    counters = {kind: SupportCounter(tree, rooted=kind in ROOTED_TREES) for kind, tree in built.items()}
    index = next(iter(counters.values())).index
    rooted = {kind: counter.rooted for kind, counter in counters.items()}
    t0 = time()
    done = 0

    def collect(job):
        nonlocal done
        for kind, found in job.result().items():
            counters[kind].add_splits(found)
        done += 1
        log.info('  %d bootstrap replicate(s) done (%s)', done, elapsed_time(time() - t0))

    with futures.ProcessPoolExecutor(max_workers=workers, mp_context=worker_context()) as executor:
        pending = deque()
        for rep_df in replicates:
            if set(rep_df.index) != set(index):
                raise ValueError('Bootstrap replicate has different samples than the reference matrix')
            pending.append(executor.submit(replicate_splits, rep_df, index, rooted, linkage, nj, me))
            # Wait when too many replicates are queued; otherwise just collect the finished ones (in order)
            while pending and (len(pending) > workers or pending[0].done()):
                collect(pending.popleft())
        while pending:
            collect(pending.popleft())

    for counter in counters.values():
        counter.assign()


def analyze_matrix(df, out_dir, name, linkage='average', nj=False, me=False, pcoa=False,
                   metadata=None, color_by=None, replicates=None, workers=1):
    """
    Build trees and the PCoA plot from a validated square distance matrix.

    :param df: square distance matrix (pandas DataFrame)
    :param out_dir: folder for the result files
    :param name: prefix of the output files
    :param replicates: optional iterable of bootstrap replicate matrices
    :param workers: number of processes building the bootstrap replicate trees
    :return: list of output files
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = list()

    built = build_trees(df, linkage, nj, me)
    if replicates is not None:
        with step('Bootstrapping with {} tree-building process(es)'.format(workers)):
            add_bootstrap_support(built, replicates, linkage, nj, me, workers)
    for kind, tree in built.items():
        out = out_dir / '{}_{}.nwk'.format(name, kind)
        trees.write_newick(tree, out)
        outputs.append(out)

    if pcoa:
        with step('Running PCoA'):
            coords, explained = ordination.pcoa(df)
            coords_file = out_dir / '{}_PCoA.tsv'.format(name)
            coords.to_csv(coords_file, sep='\t', index_label='sample', float_format='%.6g')
            meta = ordination.read_metadata(metadata, color_by) if metadata else None
            if meta is not None:
                missing = set(df.index) - set(meta.index)
                if missing:
                    log.warning('%d sample(s) are missing from the metadata file (e.g. %s)',
                                len(missing), ', '.join(sorted(missing)[:5]))
            html_file = out_dir / '{}_PCoA.html'.format(name)
            ordination.plot_pcoa(coords, explained, html_file, meta, color_by, title='PCoA of {}'.format(name))
            outputs += [coords_file, html_file]

    return outputs


class GenomeComparator:
    """The main pipeline: sample discovery, sketching, distances, then trees and PCoA."""

    def __init__(self, input_dir, output_dir, threads=1, kmer_size=21, sketch_size=10000, min_copies=2,
                 linkage='average', nj=False, me=False, pcoa=False, metadata=None, color_by=None,
                 phylip=False, force=False, clean=False, bootstrap=0):
        self.input_dir = Path(input_dir).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.sketch_dir = self.output_dir / 'sketches'
        self.tree_dir = self.output_dir / 'tree'
        self.threads = threads
        self.kmer_size = kmer_size
        self.sketch_size = sketch_size
        self.min_copies = min_copies
        self.tree_options = dict(linkage=linkage, nj=nj, me=me, pcoa=pcoa, metadata=metadata, color_by=color_by)
        self.phylip = phylip
        self.force = force
        self.clean = clean
        self.bootstrap = bootstrap

    def run(self):
        start_time = time()

        version, mash_path = mash.check_mash()
        log.info('Using Mash %s (%s), k-mer size %d, sketch size %d, %d thread(s)',
                 version, mash_path, self.kmer_size, self.sketch_size, self.threads)

        if not self.input_dir.is_dir():
            raise SampleError('Input folder "{}" does not exist or is not a folder'.format(self.input_dir))
        samples = find_samples(self.input_dir, exclude=[self.output_dir])
        n_fastq = sum(s.is_fastq for s in samples.values())
        log.info('Found %d samples (%d assemblies, %d read sets)', len(samples), len(samples) - n_fastq, n_fastq)
        if len(samples) < MIN_SAMPLES:
            raise SampleError('At least {} samples are required to build a tree, found {}'.format(
                MIN_SAMPLES, len(samples)))

        for folder in (self.output_dir, self.sketch_dir, self.tree_dir):
            folder.mkdir(parents=True, exist_ok=True)

        with step('Sketching {} samples'.format(len(samples))):
            sketches, stats, failed = self.sketch_all(samples)
        if failed:
            log.warning('%d sample(s) could not be sketched and were excluded: %s',
                        len(failed), ', '.join(sorted(failed)))
        if len(sketches) < MIN_SAMPLES:
            raise SampleError('Only {} sample(s) could be sketched, at least {} are required'.format(
                len(sketches), MIN_SAMPLES))

        all_msh = self.output_dir / 'all.msh'
        with step('Pasting sketches together'):
            mash.paste([sketches[name] for name in sorted(sketches)], all_msh)

        self.write_stats(samples, stats, failed, mash.info(all_msh), self.output_dir / 'sample_stats.tsv')

        with step('Measuring pairwise distances'):
            names, values, max_pvalue = mash.triangle(all_msh, self.threads)
            df = matrix.validate(pd.DataFrame(values, index=names, columns=names))
        if max_pvalue is not None:
            log.info('Largest Mash p-value: %g', max_pvalue)
            if max_pvalue > 0.01:
                log.warning('Some distances are not significant (p-value up to %g). '
                            'Samples may be unrelated or the sketch size too small.', max_pvalue)

        matrix_file = self.output_dir / 'all_dist.tsv'
        matrix.write_tsv(df, matrix_file)
        if self.phylip:
            matrix.write_phylip(df, self.output_dir / 'all_dist.phylip')

        replicates, workers = None, 1
        if self.bootstrap:
            sketched = [samples[name] for name in sorted(sketches)]
            replicates = self.bootstrap_matrices(sketched, self.bootstrap)
            # Mash already uses all threads; the tree-building processes share the CPUs with it.
            # Capped because each process holds a full distance matrix in memory.
            workers = max(1, min(MAX_TREE_WORKERS, self.threads // 2, self.bootstrap))
        outputs = analyze_matrix(df, self.tree_dir, 'all_dist', replicates=replicates, workers=workers,
                                 **self.tree_options)

        if self.clean:
            self.cleanup(sketches)

        log.info('Distance matrix: %s', matrix_file)
        for out in outputs:
            log.info('Output: %s', out)
        log.info('Done in %s', elapsed_time(time() - start_time))

    def sketch_all(self, samples):
        """Sketch all samples in parallel, one Mash process per sample."""
        sketches, stats, failed = dict(), dict(), dict()
        reused = 0
        with thread_pool(self.threads) as executor:
            jobs = {executor.submit(mash.sketch_or_reuse, sample, self.sketch_dir, self.kmer_size,
                                    self.sketch_size, self.min_copies, self.force): name
                    for name, sample in samples.items()}
            for done, job in enumerate(futures.as_completed(jobs), 1):
                name = jobs[job]
                try:
                    sketches[name], stats[name], was_reused = job.result()
                    reused += was_reused
                except mash.MashError as e:
                    failed[name] = str(e)
                    log.error('Could not sketch "%s": %s', name, e)
                if done % max(1, len(jobs) // 10) == 0 or done == len(jobs):
                    log.info('  %d/%d samples sketched', done, len(jobs))
        if reused:
            log.info('Reused %d existing sketch(es). Use --force to sketch everything again.', reused)
        return sketches, stats, failed

    def bootstrap_matrices(self, samples, n_replicates):
        """
        Yield one distance matrix per bootstrap replicate. Each replicate sketches the samples again with a
        different hash seed, i.e. a different random subset of k-mers. Replicate sketches are temporary.
        """
        for i in range(1, n_replicates + 1):
            t0 = time()
            seed = mash.DEFAULT_SEED + i
            with tempfile.TemporaryDirectory(dir=self.output_dir, prefix='.bootstrap_') as tmp:
                def sketch(sample):
                    prefix = Path(tmp, sample.name)
                    mash.sketch(sample, prefix, self.kmer_size, self.sketch_size, self.min_copies, seed=seed)
                    return prefix.with_name(prefix.name + '.msh')

                with thread_pool(self.threads) as executor:
                    sketch_files = list(executor.map(sketch, samples))
                all_msh = Path(tmp, 'all.msh')
                mash.paste(sketch_files, all_msh)
                names, values, _ = mash.triangle(all_msh, self.threads)
            log.debug('Bootstrap replicate %d/%d distances computed in %s', i, n_replicates, elapsed_time(time() - t0))
            yield matrix.validate(pd.DataFrame(values, index=names, columns=names))

    @staticmethod
    def write_stats(samples, sketch_stats, failed, sketch_info, stats_file):
        """
        Per-sample summary table. For assemblies "length" is the assembly size and "sequences" the number of
        contigs. For reads, "length" is the genome size estimated by Mash and "sequences" the number of reads.
        """
        with open(stats_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t', lineterminator='\n')
            writer.writerow(['sample', 'type', 'files', 'length', 'sequences', 'est_coverage', 'status'])
            for name, sample in samples.items():
                info = sketch_info.get(name, {})
                cov = sketch_stats.get(name, {}).get('est_coverage')
                writer.writerow([name, sample.file_type, len(sample.files),
                                 info.get('length', ''), info.get('num_seqs') or '',
                                 '' if cov is None else cov,
                                 'failed' if name in failed else 'ok'])

    def cleanup(self, sketches):
        """Remove only the files this tool created."""
        for msh in sketches.values():
            for f in (msh, msh.with_suffix('.json')):
                if f.exists():
                    f.unlink()
        try:
            self.sketch_dir.rmdir()
        except OSError:
            pass  # Not empty: something else lives there, leave it alone

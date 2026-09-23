"""The genome comparison pipeline: sketch -> paste -> pairwise distances -> trees / PCoA."""

import csv
import logging
from concurrent import futures
from contextlib import contextmanager
from pathlib import Path
from time import time

import pandas as pd

from . import mash, matrix, ordination, trees
from .samples import SampleError, find_samples

log = logging.getLogger(__name__)

MIN_SAMPLES = 3  # Minimum number of samples to build a tree


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


def analyze_matrix(df, out_dir, name, linkage='average', nj=False, me=False, pcoa=False,
                   metadata=None, color_by=None):
    """
    Build trees and the PCoA plot from a validated square distance matrix.

    :param df: square distance matrix (pandas DataFrame)
    :param out_dir: folder for the result files
    :param name: prefix of the output files
    :return: list of output files
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = list()

    with step('Building hierarchical clustering tree ({})'.format(linkage)):
        out = out_dir / '{}_hc.nwk'.format(name)
        trees.write_newick(trees.hc_tree(df, linkage), out)
        outputs.append(out)

    if me:
        with step('Building minimum evolution tree'):
            out = out_dir / '{}_me.nwk'.format(name)
            trees.write_newick(trees.me_tree(df), out)
            outputs.append(out)

    if nj:
        with step('Building neighbour joining tree'):
            out = out_dir / '{}_nj.nwk'.format(name)
            trees.write_newick(trees.nj_tree(df), out)
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


class MashPhylo:
    def __init__(self, input_dir, output_dir, threads=1, kmer_size=21, sketch_size=10000, min_copies=2,
                 linkage='average', nj=False, me=False, pcoa=False, metadata=None, color_by=None,
                 phylip=False, force=False, clean=False):
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

    def run(self):
        start_time = time()

        version = mash.check_mash()
        log.info('Using Mash %s, k-mer size %d, sketch size %d, %d thread(s)',
                 version, self.kmer_size, self.sketch_size, self.threads)

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

        outputs = analyze_matrix(df, self.tree_dir, 'all_dist', **self.tree_options)

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
        with futures.ThreadPoolExecutor(max_workers=self.threads) as executor:
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

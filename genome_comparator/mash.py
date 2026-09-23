"""Thin wrapper around the Mash command line tool."""

import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

MAX_KMER_SIZE = 32  # Mash limit with 64-bit hashes
DEFAULT_SEED = 42  # Mash's default hash seed
PASTE_CHUNK_SIZE = 5000  # Avoid huge single "mash paste" calls


class MashError(Exception):
    pass


def run(cmd, **kwargs):
    """Run a command, raise MashError with its stderr if it fails."""
    log.debug('Running: %s', ' '.join(map(str, cmd)))
    proc = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kwargs)
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout).strip().splitlines()
        raise MashError('Command failed (exit code {}): {}\n{}'.format(
            proc.returncode, ' '.join(map(str, cmd)), '\n'.join(detail[-10:])))
    return proc


def check_mash():
    """Make sure Mash is installed and return its version."""
    if shutil.which('mash') is None:
        raise MashError('"mash" was not found in your PATH. '
                        'Install it with "conda install -c bioconda mash".')
    version = run(['mash', '--version']).stdout.strip()
    try:
        major = int(version.split('.')[0])
    except ValueError:
        major = 0
    if major < 2:
        raise MashError('Mash version 2 or later is required (found "{}")'.format(version))
    return version


def sketch(sample, out_prefix, kmer_size, sketch_size, min_copies, seed=DEFAULT_SEED):
    """
    Sketch one sample. All files of a fastq sample (e.g. R1 and R2) are combined in a single sketch.
    The sketch ID is set to the sample name so no renaming is needed downstream.
    A different hash seed selects a different random subset of k-mers (used for bootstrapping).

    :return: dict of statistics reported by Mash (estimated coverage for reads)
    """
    cmd = ['mash', 'sketch',
           '-k', kmer_size,
           '-s', sketch_size,
           '-p', 1,
           '-S', seed,
           '-I', sample.name,
           '-o', out_prefix]
    if sample.is_fastq:
        cmd += ['-r', '-m', min_copies]
    proc = run(cmd + sample.files)

    stats = dict()
    if sample.is_fastq:
        # Mash reports "Estimated genome size: X" and "Estimated coverage: Y" on stderr for reads
        cov = re.search(r'Estimated coverage:\s*(\S+)', proc.stderr)
        if cov:
            stats['est_coverage'] = float(cov.group(1))
    return stats


def sketch_or_reuse(sample, sketch_dir, kmer_size, sketch_size, min_copies, force=False):
    """
    Sketch a sample unless an up-to-date sketch made with the same parameters already exists.
    A small json file is saved next to each sketch to remember how it was made.

    :return: (sketch_path, stats, reused)
    """
    prefix = Path(sketch_dir, sample.name)
    msh = prefix.with_name(prefix.name + '.msh')
    meta_file = prefix.with_name(prefix.name + '.json')
    params = {'kmer_size': kmer_size, 'sketch_size': sketch_size,
              'min_copies': min_copies if sample.is_fastq else None,
              'files': [os.path.abspath(f) for f in sample.files]}

    if not force and msh.exists() and meta_file.exists():
        try:
            meta = json.loads(meta_file.read_text())
            newest_input = max(os.path.getmtime(f) for f in sample.files)
            if meta.get('params') == params and msh.stat().st_mtime >= newest_input:
                return msh, meta.get('stats', {}), True
        except (ValueError, OSError):
            pass  # Corrupted metadata, just sketch again

    # Make sure no stale sketch survives if sketching fails
    for f in (msh, meta_file):
        if f.exists():
            f.unlink()
    stats = sketch(sample, prefix, kmer_size, sketch_size, min_copies)
    meta_file.write_text(json.dumps({'params': params, 'stats': stats}))
    return msh, stats, False


def paste(sketch_files, output_msh):
    """
    Combine many sketch files into a single one.
    Large inputs are pasted in chunks, then the chunks are pasted together.
    """
    output_msh = Path(output_msh)
    if output_msh.exists():
        output_msh.unlink()  # "mash paste" refuses to overwrite an existing file

    with tempfile.TemporaryDirectory(dir=output_msh.parent, prefix='.paste_') as tmp:
        files = [str(f) for f in sketch_files]
        level = 0
        while len(files) > PASTE_CHUNK_SIZE:
            chunks = [files[i:i + PASTE_CHUNK_SIZE] for i in range(0, len(files), PASTE_CHUNK_SIZE)]
            files = [_paste_list(chunk, Path(tmp, 'chunk{}_{}.msh'.format(level, i)))
                     for i, chunk in enumerate(chunks)]
            level += 1
        _paste_list(files, output_msh)


def _paste_list(files, output_msh):
    with tempfile.NamedTemporaryFile('w', dir=output_msh.parent, suffix='.list', delete=False) as fh:
        fh.write('\n'.join(files) + '\n')
    try:
        run(['mash', 'paste', output_msh, '-l', fh.name])
    finally:
        os.unlink(fh.name)
    return str(output_msh)


def info(msh):
    """
    Read per-sketch information from a sketch file.

    :return: dict {sketch_id: {'length': int, 'num_seqs': int or None}}
             'length' is the total sequence length for assemblies and the estimated genome size for reads.
    """
    result = dict()
    for line in run(['mash', 'info', '-t', msh]).stdout.splitlines():
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        if len(fields) < 3:
            continue
        seqs = re.match(r'\[(\d+) seqs?\]', fields[3]) if len(fields) > 3 else None
        result[fields[2]] = {'length': int(fields[1]),
                             'num_seqs': int(seqs.group(1)) if seqs else None}
    return result


def triangle(msh, threads=1):
    """
    Compute all pairwise distances with "mash triangle" (only half of the comparisons are needed).

    :return: (names, square symmetric numpy distance matrix, max p-value or None)
    """
    with tempfile.TemporaryFile('w+') as err:
        cmd = ['mash', 'triangle', '-p', str(threads), str(msh)]
        log.debug('Running: %s', ' '.join(cmd))
        parse_error = None
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, text=True) as proc:
            try:
                names, matrix = parse_triangle(proc.stdout)
            except MashError as e:
                parse_error = e
        err.seek(0)
        stderr = err.read()
    # Report Mash's own error first, the parsing error is usually just a consequence of it
    if proc.returncode != 0:
        raise MashError('mash triangle failed (exit code {}):\n{}'.format(proc.returncode, stderr.strip()))
    if parse_error:
        raise parse_error

    pvalue = re.search(r'Max p-value:\s*(\S+)', stderr)
    return names, matrix, float(pvalue.group(1)) if pvalue else None


def parse_triangle(lines):
    """
    Parse "mash triangle" output (streamed line by line to keep memory low on big datasets):

        <tab>3
        A
        B<tab>0.1
        C<tab>0.2<tab>0.3
    """
    lines = iter(lines)
    try:
        n = int(next(lines).strip())
    except (StopIteration, ValueError):
        raise MashError('Unexpected "mash triangle" output')

    names = list()
    matrix = np.zeros((n, n), dtype=float)
    for i, line in enumerate(lines):
        fields = line.rstrip('\n').split('\t')
        if i >= n or len(fields) != i + 1:
            raise MashError('Malformed "mash triangle" output at row {}'.format(i + 1))
        names.append(fields[0])
        if i:
            row = np.array(fields[1:], dtype=float)
            matrix[i, :i] = row
            matrix[:i, i] = row
    if len(names) != n:
        raise MashError('Truncated "mash triangle" output: {} of {} rows'.format(len(names), n))
    return names, matrix

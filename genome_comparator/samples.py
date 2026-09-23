"""Discover input sequence files and group them into samples."""

import os
import re
from dataclasses import dataclass, field
from pathlib import Path

FASTA_EXTENSIONS = ('.fna', '.fasta', '.fa', '.fas')
FASTQ_EXTENSIONS = ('.fastq', '.fq')
COMPRESSION_EXTENSIONS = ('', '.gz')

# Paired-end / multi-file suffixes stripped from fastq file names to get the sample name.
# e.g. "S1_R1.fastq.gz", "S1_S12_L001_R2_001.fastq.gz", "SRR123_1.fastq.gz"
READ_SUFFIX = re.compile(r'_R?[12](_001)?$')


class SampleError(Exception):
    pass


@dataclass
class Sample:
    name: str
    file_type: str  # 'fasta' or 'fastq'
    files: list = field(default_factory=list)

    @property
    def is_fastq(self):
        return self.file_type == 'fastq'


def split_extension(filename):
    """
    Return (stem, file_type) for a supported sequence file, or None if the file is not supported.

    Only the known extension is removed, so dots elsewhere in the name are kept
    (e.g. "E.coli_K12.v2.fasta.gz" -> "E.coli_K12.v2").
    """
    lower = filename.lower()
    for file_type, extensions in (('fasta', FASTA_EXTENSIONS), ('fastq', FASTQ_EXTENSIONS)):
        for ext in extensions:
            for comp in COMPRESSION_EXTENSIONS:
                suffix = ext + comp
                if lower.endswith(suffix) and len(filename) > len(suffix):
                    return filename[:-len(suffix)], file_type
    return None


def sample_name_from_file(filename):
    """Return (sample_name, file_type) or None if the file is not a supported sequence file."""
    parsed = split_extension(filename)
    if parsed is None:
        return None
    stem, file_type = parsed
    if file_type == 'fastq':
        stem = READ_SUFFIX.sub('', stem) or stem
    return stem, file_type


def find_samples(input_folder, exclude=()):
    """
    Recursively look for fasta/fastq files and group them per sample.

    Fastq files sharing a sample name (e.g. R1/R2) are grouped together and sketched as one sample.
    Every other name collision is an error, because silently merging different genomes would
    produce a wrong distance matrix.

    :param input_folder: folder to search recursively
    :param exclude: folders to skip (e.g. the output folder if it is nested in the input folder)
    :return: dict {sample_name: Sample}, sorted by sample name
    """
    exclude = {Path(p).resolve() for p in exclude}
    samples = dict()
    errors = list()

    for root, dirs, filenames in os.walk(input_folder):
        # Prune excluded folders and keep traversal order deterministic
        dirs[:] = sorted(d for d in dirs if Path(root, d).resolve() not in exclude)
        for filename in sorted(filenames):
            parsed = sample_name_from_file(filename)
            if parsed is None:
                continue  # Ignore other files
            name, file_type = parsed
            file_path = os.path.join(root, filename)

            sample = samples.get(name)
            if sample is None:
                samples[name] = Sample(name, file_type, [file_path])
            elif sample.file_type != file_type:
                errors.append('Sample "{}" has both fasta and fastq files: {}'.format(
                    name, ', '.join(sample.files + [file_path])))
            elif file_type == 'fasta':
                errors.append('Sample name "{}" is shared by several fasta files: {}'.format(
                    name, ', '.join(sample.files + [file_path])))
            elif filename in {os.path.basename(f) for f in sample.files}:
                errors.append('File "{}" found in several folders: {}'.format(
                    filename, ', '.join(sample.files + [file_path])))
            else:
                sample.files.append(file_path)

    if errors:
        raise SampleError('Ambiguous sample names. Rename or remove the offending files:\n  '
                          + '\n  '.join(errors))

    return dict(sorted(samples.items()))

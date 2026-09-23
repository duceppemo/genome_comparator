"""End-to-end tests. Skipped when Mash is not installed."""

import random
import shutil

import pandas as pd
import pytest

from genome_comparator.cli import main
from genome_comparator.trees import read_newick

pytestmark = pytest.mark.skipif(shutil.which('mash') is None, reason='mash is not installed')


def mutate(seq, rate, rng):
    return ''.join(rng.choice('ACGT') if rng.random() < rate else c for c in seq)


@pytest.fixture(scope='module')
def genomes(tmp_path_factory):
    """Four related synthetic assemblies and one paired-end read set."""
    rng = random.Random(1)
    folder = tmp_path_factory.mktemp('genomes')
    ref = ''.join(rng.choice('ACGT') for _ in range(50000))
    for name, rate in (('A', 0), ('B', 0.01), ('C', 0.03), ('Iso_R10', 0.05)):
        seq = mutate(ref, rate, rng)
        (folder / '{}.fasta'.format(name)).write_text('>c1\n{}\n>c2\n{}\n'.format(seq[:25000], seq[25000:]))
    reads = mutate(ref, 0.02, rng)
    for mate in ('R1', 'R2'):
        with open(folder / 'R_{}.fastq'.format(mate), 'w') as f:
            for i in range(3000):
                p = rng.randint(0, len(reads) - 100)
                f.write('@r{}\n{}\n+\n{}\n'.format(i, reads[p:p + 100], 'I' * 100))
    return folder


def test_full_pipeline(genomes, tmp_path):
    out = tmp_path / 'out'
    main(['-i', str(genomes), '-o', str(out), '-t', '2', '-s', '1000', '--nj', '--me', '--pcoa', '--phylip'])

    df = pd.read_csv(out / 'all_dist.tsv', sep='\t', index_col=0)
    assert list(df.index) == list(df.columns) == ['A', 'B', 'C', 'Iso_R10', 'R']
    assert df.loc['A', 'B'] < df.loc['A', 'C'] < df.loc['A', 'Iso_R10']

    stats = pd.read_csv(out / 'sample_stats.tsv', sep='\t', index_col=0)
    assert stats.loc['A', 'length'] == 50000 and stats.loc['A', 'sequences'] == 2
    assert stats.loc['R', 'type'] == 'fastq' and stats.loc['R', 'files'] == 2
    assert stats.loc['R', 'est_coverage'] > 0

    for suffix in ('hc.nwk', 'nj.nwk', 'me.nwk', 'PCoA.html', 'PCoA.tsv'):
        assert (out / 'tree' / 'all_dist_{}'.format(suffix)).exists()
    tree = read_newick(out / 'tree' / 'all_dist_hc.nwk')
    assert sorted(t.name for t in tree.tips()) == list(df.index)


def test_rerun_reuses_sketches_and_updates_results(genomes, tmp_path, caplog):
    out = tmp_path / 'out'
    main(['-i', str(genomes), '-o', str(out), '-t', '2', '-s', '1000'])
    caplog.clear()
    # A second run in the same folder used to silently keep the old all.msh
    main(['-i', str(genomes), '-o', str(out), '-t', '2', '-s', '1000'])
    assert 'Reused 5 existing sketch(es)' in caplog.text

    # Changing a parameter must invalidate the sketches
    caplog.clear()
    main(['-i', str(genomes), '-o', str(out), '-t', '2', '-s', '500', '--clean'])
    assert 'Reused' not in caplog.text
    assert not (out / 'sketches').exists()
    assert (out / 'all.msh').exists()


def test_bad_file_is_excluded(genomes, tmp_path, caplog):
    folder = tmp_path / 'in'
    shutil.copytree(genomes, folder)
    (folder / 'broken.fasta').write_text('not a fasta file\n')
    out = tmp_path / 'out'
    main(['-i', str(folder), '-o', str(out), '-t', '2', '-s', '1000'])
    assert 'could not be sketched' in caplog.text
    stats = pd.read_csv(out / 'sample_stats.tsv', sep='\t', index_col=0)
    assert stats.loc['broken', 'status'] == 'failed'
    assert 'broken' not in pd.read_csv(out / 'all_dist.tsv', sep='\t', index_col=0).index


def test_too_few_samples_exits_with_error(tmp_path):
    (tmp_path / 'A.fasta').write_text('>a\nACGT\n')
    with pytest.raises(SystemExit) as e:
        main(['-i', str(tmp_path), '-o', str(tmp_path / 'out')])
    assert e.value.code == 1

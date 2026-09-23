import io

import numpy as np
import pandas as pd
import pytest
from skbio import TreeNode

from genome_comparator import matrix, trees
from genome_comparator.mash import MashError, parse_triangle
from genome_comparator.pipeline import elapsed_time
from genome_comparator.samples import SampleError, find_samples, sample_name_from_file
from genome_comparator.tree_tools import collapse, rename_tips


@pytest.mark.parametrize('filename, expected', [
    ('S1.fasta', ('S1', 'fasta')),
    ('E.coli_K12.v2.fna.gz', ('E.coli_K12.v2', 'fasta')),
    ('Strain_R1_final.fasta', ('Strain_R1_final', 'fasta')),  # Fasta names are never trimmed
    ('Iso_R10.fasta', ('Iso_R10', 'fasta')),
    ('S1_R1.fastq.gz', ('S1', 'fastq')),
    ('S1_S12_L001_R2_001.fastq.gz', ('S1_S12_L001', 'fastq')),
    ('SRR123_1.fq', ('SRR123', 'fastq')),
    ('Iso_R10.fastq', ('Iso_R10', 'fastq')),  # "_R10" is not a read suffix
    ('notes.txt', None),
    ('.fasta', None),
])
def test_sample_name_from_file(filename, expected):
    assert sample_name_from_file(filename) == expected


def touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('')


def test_find_samples_groups_reads(tmp_path):
    for f in ('A_R1.fastq.gz', 'A_R2.fastq.gz', 'B.fasta', 'sub/C.fa', 'readme.md'):
        touch(tmp_path / f)
    samples = find_samples(tmp_path)
    assert list(samples) == ['A', 'B', 'C']
    assert samples['A'].is_fastq and len(samples['A'].files) == 2


def test_find_samples_excludes_output(tmp_path):
    touch(tmp_path / 'A.fasta')
    touch(tmp_path / 'out' / 'B.fasta')
    assert list(find_samples(tmp_path, exclude=[tmp_path / 'out'])) == ['A']


@pytest.mark.parametrize('files', [
    ('A.fasta', 'A.fna'),  # Two assemblies with the same name
    ('A.fasta', 'A_R1.fastq'),  # Mixed types
    ('x/A_R1.fastq', 'y/A_R1.fastq'),  # Same file in two folders
])
def test_find_samples_rejects_ambiguous_names(tmp_path, files):
    for f in files:
        touch(tmp_path / f)
    with pytest.raises(SampleError):
        find_samples(tmp_path)


def test_parse_triangle():
    names, m = parse_triangle(['\t3\n', 'A\n', 'B\t0.1\n', 'C\t0.2\t0.3\n'])
    assert names == ['A', 'B', 'C']
    np.testing.assert_allclose(m, [[0, .1, .2], [.1, 0, .3], [.2, .3, 0]])


def test_parse_triangle_truncated():
    with pytest.raises(MashError):
        parse_triangle(['\t3\n', 'A\n', 'B\t0.1\n'])


def square(names, values):
    return pd.DataFrame(np.array(values, dtype=float), index=names, columns=names)


def test_validate_sorts_rows_and_columns():
    df = pd.DataFrame([[0, .2, .1], [.2, 0, .3], [.1, .3, 0]], index=['B', 'C', 'A'], columns=['B', 'C', 'A'])
    out = matrix.validate(df)
    assert list(out.index) == list(out.columns) == ['A', 'B', 'C']
    assert out.loc['A', 'C'] == pytest.approx(.3)


@pytest.mark.parametrize('df', [
    square(['A', 'B'], [[0, .1], [.2, 0]]),  # Asymmetric
    square(['A', 'A'], [[0, .1], [.1, 0]]),  # Duplicated names
    pd.DataFrame([[0, .1]], index=['A'], columns=['A', 'B']),  # Not square
    square(['A', 'B'], [[0, np.nan], [np.nan, 0]]),  # Missing values
    square(['A', 'B'], [[0, np.inf], [np.inf, 0]]),  # Infinite values
    square(['A', 'B'], [[0, -.1], [-.1, 0]]),  # Negative values
])
def test_validate_rejects_bad_matrices(df):
    with pytest.raises(matrix.MatrixError):
        matrix.validate(df)


def test_read_matrix_keeps_numeric_looking_names(tmp_path):
    # "001" used to be parsed as the number 1, which no longer matched the "001" column name
    path = tmp_path / 'm.tsv'
    path.write_text('#query\t001\t002\t1e3\n001\t0\t0.1\t0.2\n002\t0.1\t0\t0.3\n1e3\t0.2\t0.3\t0\n')
    assert list(matrix.read_matrix(path).index) == ['001', '002', '1e3']
    path = tmp_path / 'm.csv'
    path.write_text(',7,8,9\n7,0,0.1,0.2\n8,0.1,0,0.3\n9,0.2,0.3,0\n')
    assert list(matrix.read_matrix(path).index) == ['7', '8', '9']


def test_validate_writes_integral_names_without_decimals():
    # Excel hands numeric sample IDs over as floats
    df = pd.DataFrame(np.zeros((2, 2)), index=[1.0, 2.0], columns=[1.0, 2.0])
    assert list(matrix.validate(df).index) == ['1', '2']


def test_read_write_matrix_roundtrip(tmp_path):
    df = square(['A', 'B', 'C'], [[0, .01, .02], [.01, 0, .03], [.02, .03, 0]])
    matrix.write_tsv(df, tmp_path / 'm.tsv')
    assert (tmp_path / 'm.tsv').read_text().startswith('#query\tA\tB\tC\n')
    pd.testing.assert_frame_equal(matrix.read_matrix(tmp_path / 'm.tsv'), df)


def test_hc_tree_keeps_small_branch_lengths():
    # The old implementation rounded branch lengths to 2 decimals, turning Mash distances into 0.00
    df = square(['A', 'B', 'C'], [[0, .002, .01], [.002, 0, .01], [.01, .01, 0]])
    tree = trees.hc_tree(df, 'average')
    assert tree.find('A').length == pytest.approx(.001)
    assert tree.find('A').distance(tree.find('B')) == pytest.approx(.002)


def test_newick_quotes_all_labels():
    tree = TreeNode.read(io.StringIO("(('it''s x':0.1,'a,b':0.2):0.3,c_d:0.4);"), convert_underscores=False)
    nwk = trees.to_newick(tree)
    assert nwk == "(('it''s x':0.1,'a,b':0.2):0.3,'c_d':0.4);\n"
    back = TreeNode.read(io.StringIO(nwk), convert_underscores=False)
    assert sorted(t.name for t in back.tips()) == sorted(["it's x", 'a,b', 'c_d'])


def test_newick_deep_tree_has_no_recursion_limit():
    nwk = '(' * 5000 + 'A' + ''.join(',B{}:1)'.format(i) for i in range(5000)) + ';'
    tree = TreeNode.read(io.StringIO(nwk))
    assert trees.to_newick(tree).count('(') == 5000


def test_collapse():
    tree = TreeNode.read(io.StringIO('((A:0.001,B:0.001):0.1,(C:0.2,D:0.2):0.1);'))
    assert collapse(tree, 0.01) == 1
    assert sorted(t.name for t in tree.tips()) == ['A {B}', 'C', 'D']


def test_rename_exact_match_only():
    tree = TreeNode.read(io.StringIO('(S1:1,S10:1,S11:1);'))
    missing, duplicates = rename_tips(tree, {'S1': 'new', 'S2': 'x'})
    assert sorted(t.name for t in tree.tips()) == ['S10', 'S11', 'new']
    assert missing == {'S2'} and not duplicates


def test_rename_reports_duplicate_names():
    tree = TreeNode.read(io.StringIO('(S1:1,S2:1,S3:1);'))
    _, duplicates = rename_tips(tree, {'S1': 'same', 'S2': 'same'})
    assert duplicates == {'same'}


@pytest.mark.parametrize('seconds, expected', [(0, '0s'), (59.6, '1m'), (3725, '1h2m5s')])
def test_elapsed_time(seconds, expected):
    assert elapsed_time(seconds) == expected


def read(nwk):
    return TreeNode.read(io.StringIO(nwk), convert_underscores=False)


def test_rooted_clades_differ_from_unrooted_splits():
    from genome_comparator.bootstrap import node_splits
    index = {n: i for i, n in enumerate('ABCD')}
    # Same unrooted tree AB|CD, rooted differently
    t1, t2 = read('((A,B),(C,D));'), read('(A,(B,(C,D)));')
    assert set(node_splits(t1, index, rooted=False).values()) == set(node_splits(t2, index, rooted=False).values())
    assert set(node_splits(t1, index, rooted=True).values()) != set(node_splits(t2, index, rooted=True).values())


def test_support_counter():
    from genome_comparator.bootstrap import SupportCounter
    ref = read('(((A,B),C),(D,E));')
    counter = SupportCounter(ref, rooted=True)
    counter.add(read('(((A,B),C),(D,E));'))
    counter.add(read('(((A,C),B),(D,E));'))
    counter.assign()
    supports = {frozenset(t.name for t in n.tips()): n.support for n in ref.non_tips()}
    assert supports == {frozenset('AB'): 50, frozenset('ABC'): 100, frozenset('DE'): 100}
    assert trees.to_newick(ref) == "((('A','B')50,'C')100,('D','E')100);\n"


def test_bootstrap_support_in_parallel():
    from genome_comparator.pipeline import add_bootstrap_support, build_trees
    names = list('ABCDE')
    ref = square(names, [[0, 1, 4, 8, 8], [1, 0, 4, 8, 8], [4, 4, 0, 8, 8], [8, 8, 8, 0, 2], [8, 8, 8, 2, 0]])
    swapped = ref.copy()  # A is now closer to C than to B
    swapped.loc['A', 'C'] = swapped.loc['C', 'A'] = 0.5
    built = build_trees(ref, nj=True, me=True, verbose=False)
    add_bootstrap_support(built, [ref, ref, ref, swapped], nj=True, me=True, workers=2)
    hc = {frozenset(t.name for t in n.tips()): n.support for n in built['hc'].non_tips()}
    assert hc[frozenset('AB')] == 75 and hc[frozenset('DE')] == 100
    for kind in ('nj', 'me'):
        supports = [n.support for n in built[kind].non_tips()]
        assert supports and all(0 <= s <= 100 for s in supports)


def test_thread_pool_drops_queued_jobs_on_interrupt():
    import time
    from genome_comparator.pipeline import thread_pool
    ran = list()

    def job(i):
        time.sleep(0.01)
        ran.append(i)

    with pytest.raises(KeyboardInterrupt):
        with thread_pool(1) as executor:
            [executor.submit(job, i) for i in range(50)]
            raise KeyboardInterrupt
    assert len(ran) < 50  # Without cancellation, every queued job runs before the executor shuts down


@pytest.fixture
def fake_mash(tmp_path, monkeypatch):
    """A fake "mash" next to a fake Python interpreter, and an empty PATH."""
    from genome_comparator import mash
    env_bin = tmp_path / 'env' / 'bin'
    env_bin.mkdir(parents=True)
    fake = env_bin / 'mash'
    fake.write_text('#!/bin/sh\necho 2.3\n')
    fake.chmod(0o755)
    monkeypatch.setenv('PATH', str(tmp_path / 'empty'))
    monkeypatch.setattr(mash.sys, 'executable', str(env_bin / 'python'))
    mash.executable.cache_clear()
    yield fake
    mash.executable.cache_clear()


def test_mash_fallback_to_env_bin(fake_mash):
    from genome_comparator import mash
    assert mash.executable() == str(fake_mash)
    assert mash.check_mash() == ('2.3', str(fake_mash))
    assert mash.command(['mash', 'info', 1]) == [str(fake_mash), 'info', '1']


def test_mash_in_path_wins(fake_mash, tmp_path, monkeypatch):
    from genome_comparator import mash
    path_bin = tmp_path / 'path_bin'
    path_bin.mkdir()
    in_path = path_bin / 'mash'
    in_path.write_text('#!/bin/sh\necho 2.2\n')
    in_path.chmod(0o755)
    monkeypatch.setenv('PATH', str(path_bin))
    assert mash.executable() == str(in_path)


def test_mash_not_found(fake_mash):
    from genome_comparator import mash
    fake_mash.unlink()
    assert mash.executable() is None
    with pytest.raises(MashError, match='not found'):
        mash.check_mash()


def test_category_styles_are_unique_and_stable():
    from genome_comparator.ordination import NEUTRAL, category_styles
    values = pd.Series(['b', 'a', None, 'c', 'a'])
    labels, styles, order = category_styles(values)
    assert order == ['a', 'b', 'c', 'Unknown']
    assert list(labels) == ['b', 'a', 'Unknown', 'c', 'a']
    assert styles['Unknown'][0] == NEUTRAL
    assert len({styles[c] for c in 'abc'}) == 3
    # Adding a category does not change the style of the others
    _, styles2, _ = category_styles(pd.Series(['a', 'b', 'c', 'aa']))
    assert styles2['a'] == styles['a']


def test_category_styles_fold_extra_categories_into_other():
    from genome_comparator.ordination import PALETTE, SYMBOLS, category_styles
    n = len(PALETTE) * len(SYMBOLS) + 5
    labels, styles, order = category_styles(pd.Series(['c{:02d}'.format(i) for i in range(n)] + ['c00'] * 3))
    assert order[-1] == 'Other' and 'c00' in order  # The most frequent category is kept
    real = [c for c in order if c != 'Other']
    assert len({styles[c] for c in real}) == len(real) == len(PALETTE) * len(SYMBOLS) - 1


def test_category_named_other_keeps_its_own_style():
    from genome_comparator.ordination import NEUTRAL, category_styles
    _, styles, order = category_styles(pd.Series(['x', 'Other', 'y']))
    assert order == ['Other', 'x', 'y']  # Listed once, as a normal category
    assert styles['Other'][0] != NEUTRAL


def test_metadata_rejects_duplicated_sample_names(tmp_path):
    from genome_comparator.ordination import read_metadata
    path = tmp_path / 'meta.tsv'
    path.write_text('sample\tgroup\nA\tx\nA\ty\nB\tz\n')
    with pytest.raises(ValueError, match='A'):
        read_metadata(path, 'group')


def test_pcoa_html(tmp_path):
    from genome_comparator import ordination
    df = square(list('ABCD'), [[0, .1, .5, .5], [.1, 0, .5, .5], [.5, .5, 0, .1], [.5, .5, .1, 0]])
    coords, explained = ordination.pcoa(df)
    meta = pd.DataFrame({'group': ['x', 'x', 'y']}, index=['A', 'B', 'C'])  # D missing
    fig = ordination.pcoa_figure(coords, explained, meta, 'group')
    assert {t.name for t in fig.data} == {'x', 'y', 'Unknown'}
    ordination.plot_pcoa(coords, explained, tmp_path / 'p.html', meta, 'group')
    assert (tmp_path / 'p.html').stat().st_size > 0


def test_parse_info_counts_single_sequence():
    from genome_comparator.mash import parse_info
    text = ('#Hashes\tLength\tID\tComment\n'
            '1000\t2905187\tF2365\tNC_002973.6 Listeria monocytogenes F2365, complete sequence\n'
            '1000\t3091600\tR2-502\t[2 seqs] NC_021838.1 Listeria monocytogenes [...]\n'
            '1000\t191866\treads\t[8000 seqs] r0 [...]\n')
    assert parse_info(text) == {'F2365': {'length': 2905187, 'num_seqs': 1},
                                'R2-502': {'length': 3091600, 'num_seqs': 2},
                                'reads': {'length': 191866, 'num_seqs': 8000}}


def test_mash_phylo_is_a_deprecated_alias(capsys):
    from genome_comparator.cli import mash_phylo_main
    with pytest.raises(SystemExit) as e:
        mash_phylo_main(['--version'])
    assert e.value.code == 0
    out = capsys.readouterr()
    assert 'deprecated' in out.err and 'genome-comparator' in out.err
    assert out.out.startswith('genome-comparator ')


def test_python_m_runs_main_command():
    import subprocess
    import sys
    proc = subprocess.run([sys.executable, '-m', 'genome_comparator', '--version'], capture_output=True, text=True)
    assert proc.returncode == 0 and proc.stdout.startswith('genome-comparator ')


def test_me_tree_with_three_samples():
    df = square(list('ABC'), [[0, .01, .05], [.01, 0, .06], [.05, .06, 0]])
    assert sorted(t.name for t in trees.me_tree(df).tips()) == ['A', 'B', 'C']


def test_support_values_survive_read_write(tmp_path):
    path = tmp_path / 't.nwk'
    path.write_text("(('A':0.1,'B':0.1)95:0.2,'C':0.3,'D':0.1);\n")
    tree = trees.read_newick(path)
    assert trees.to_newick(tree) == "(('A':0.1,'B':0.1)95:0.2,'C':0.3,'D':0.1);\n"
    collapse(tree, 0.5)  # Collapsed clades become tips named after their tips, not after their support
    assert "'A {B}'" in trees.to_newick(tree)


def test_pcoa_of_identical_genomes():
    from genome_comparator import ordination
    _, explained = ordination.pcoa(square(list('ABCD'), np.zeros((4, 4))))
    assert explained == [0.0, 0.0, 0.0]


def test_xls_without_xlrd(tmp_path, monkeypatch):
    def missing(*args, **kwargs):
        raise ImportError('Missing optional dependency xlrd')
    monkeypatch.setattr(pd, 'read_excel', missing)
    with pytest.raises(matrix.MatrixError, match='xlrd'):
        matrix.read_matrix(tmp_path / 'm.xls')

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
])
def test_validate_rejects_bad_matrices(df):
    with pytest.raises(matrix.MatrixError):
        matrix.validate(df)


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
    missing = rename_tips(tree, {'S1': 'new', 'S2': 'x'})
    assert sorted(t.name for t in tree.tips()) == ['S10', 'S11', 'new']
    assert missing == {'S2'}


@pytest.mark.parametrize('seconds, expected', [(0, '0s'), (59.6, '1m'), (3725, '1h2m5s')])
def test_elapsed_time(seconds, expected):
    assert elapsed_time(seconds) == expected

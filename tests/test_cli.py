"""Command line entry points: argument checks, error handling and the tools that work without Mash."""

import logging
from argparse import ArgumentTypeError

import pytest

from genome_comparator import cli, trees
from genome_comparator.cli import collapse_main, dendrogram_main, rename_main
from genome_comparator.matrix import MatrixError
from genome_comparator.tree_tools import read_rename_table

MATRIX = ('\tA\tB\tC\tD\n'
          'A\t0\t0.01\t0.05\t0.06\n'
          'B\t0.01\t0\t0.05\t0.06\n'
          'C\t0.05\t0.05\t0\t0.02\n'
          'D\t0.06\t0.06\t0.02\t0\n')


@pytest.fixture
def matrix_file(tmp_path):
    path = tmp_path / 'matrix.tsv'
    path.write_text(MATRIX)
    return path


@pytest.fixture
def tree_file(tmp_path):
    path = tmp_path / 'tree.nwk'
    path.write_text("(('A':0.001,'B':0.001)95:0.1,('C':0.02,'D':0.02)80:0.1);\n")
    return path


def exit_code(func, argv):
    with pytest.raises(SystemExit) as e:
        func(argv)
    return e.value.code


def tip_names(path):
    return sorted(t.name for t in trees.read_newick(path).tips())


@pytest.mark.parametrize('check, value, message', [
    (cli.int_range(1), 'x', 'not an integer'),
    (cli.int_range(1), '0', 'must be >= 1'),
    (cli.int_range(1, 32), '33', 'between 1 and 32'),
])
def test_int_range_rejects(check, value, message):
    with pytest.raises(ArgumentTypeError, match=message):
        check(value)


def test_int_range_accepts():
    assert cli.int_range(1, 32)('32') == 32


def test_available_cpus_without_sched_getaffinity(monkeypatch):
    monkeypatch.delattr(cli.os, 'sched_getaffinity', raising=False)
    monkeypatch.setattr(cli.os, 'cpu_count', lambda: None)
    assert cli.available_cpus() == 1


def test_threads_capped_to_available_cpus(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, 'available_cpus', lambda: 2)
    seen = dict()

    class FakeComparator:
        def __init__(self, *args, **kwargs):
            seen.update(kwargs)

        def run(self):
            pass

    monkeypatch.setattr(cli, 'GenomeComparator', FakeComparator)
    cli.main(['-i', str(tmp_path), '-o', str(tmp_path / 'out'), '-t', '64'])
    assert seen['threads'] == 2
    assert (tmp_path / 'out' / 'genome_comparator.log').exists()


@pytest.mark.parametrize('extra, message', [
    (['--pcoa', '--color-by', 'species'], '--color-by requires --metadata'),
    (['--metadata', 'meta.tsv'], '--metadata is only used with --pcoa'),
    (['--pcoa', '--metadata', 'missing.tsv'], 'missing.tsv'),
])
def test_tree_arguments_checked_before_running(matrix_file, tmp_path, capsys, extra, message):
    assert exit_code(dendrogram_main, ['-i', str(matrix_file), '-o', str(tmp_path / 'out'), *extra]) == 2
    assert message in capsys.readouterr().err
    assert not (tmp_path / 'out').exists()


def test_color_by_unknown_column(matrix_file, tmp_path, capsys):
    meta = tmp_path / 'meta.tsv'
    meta.write_text('sample\tspecies\nA\tx\n')
    argv = ['-i', str(matrix_file), '-o', str(tmp_path / 'out'), '--pcoa', '--metadata', str(meta),
            '--color-by', 'serotype']
    assert exit_code(dendrogram_main, argv) == 2
    assert 'Column "serotype" not found' in capsys.readouterr().err


def test_run_safely_reports_user_errors(caplog):
    def fail():
        raise MatrixError('bad matrix')
    with pytest.raises(SystemExit) as e:
        cli.run_safely(fail)
    assert e.value.code == 1
    assert 'bad matrix' in caplog.text


def test_run_safely_interrupted(caplog):
    def interrupt():
        raise KeyboardInterrupt
    with pytest.raises(SystemExit) as e:
        cli.run_safely(interrupt)
    assert e.value.code == 130
    assert 'Interrupted' in caplog.text


def test_setup_logging_replaces_handlers(tmp_path):
    cli.setup_logging(log_file=tmp_path / 'first.log')
    cli.setup_logging(verbose=True)
    log = logging.getLogger('genome_comparator')
    assert len(log.handlers) == 1 and log.level == logging.DEBUG


def test_dendrogram_from_matrix(matrix_file, tmp_path):
    meta = tmp_path / 'meta.tsv'
    meta.write_text('sample\tspecies\nA\tx\nB\tx\nC\ty\n')
    out = tmp_path / 'out'
    dendrogram_main(['-i', str(matrix_file), '-o', str(out), '--nj', '--me', '--pcoa',
                     '--metadata', str(meta), '--color-by', 'species'])
    for kind in ('hc', 'nj', 'me'):
        assert tip_names(out / 'matrix_{}.nwk'.format(kind)) == ['A', 'B', 'C', 'D']
    assert (out / 'matrix_PCoA.tsv').exists()
    assert 'species' in (out / 'matrix_PCoA.html').read_text()


def test_dendrogram_needs_enough_samples(tmp_path, caplog):
    path = tmp_path / 'small.tsv'
    path.write_text('\tA\tB\nA\t0\t0.1\nB\t0.1\t0\n')
    assert exit_code(dendrogram_main, ['-i', str(path), '-o', str(tmp_path / 'out')]) == 1
    assert 'At least' in caplog.text


def test_dendrogram_rejects_bad_matrix(tmp_path, caplog):
    path = tmp_path / 'bad.tsv'
    path.write_text('\tA\tB\tC\nA\t0\t0.1\t0.2\nB\t0.1\t0\t0.3\n')
    assert exit_code(dendrogram_main, ['-i', str(path), '-o', str(tmp_path / 'out')]) == 1
    assert 'not square' in caplog.text


def test_dendrogram_missing_input(tmp_path, caplog):
    assert exit_code(dendrogram_main, ['-i', str(tmp_path / 'none.tsv'), '-o', str(tmp_path / 'out')]) == 1
    assert 'none.tsv' in caplog.text


def test_tree_collapser(tree_file, tmp_path, caplog):
    out = tmp_path / 'collapsed.nwk'
    with caplog.at_level(logging.INFO, logger='genome_comparator'):
        collapse_main(['-i', str(tree_file), '-o', str(out), '-d', '0.01'])
    assert tip_names(out) == ['A {B}', 'C', 'D']
    assert '80' in out.read_text()  # Support of the remaining clade is kept
    assert 'Collapsed 1 clade(s)' in caplog.text


def test_tree_collapser_missing_input(tmp_path, caplog):
    argv = ['-i', str(tmp_path / 'none.nwk'), '-o', str(tmp_path / 'out.nwk'), '-d', '0.01']
    assert exit_code(collapse_main, argv) == 1
    assert 'none.nwk' in caplog.text


def test_tree_renamer(tree_file, tmp_path, caplog):
    table = tmp_path / 'rename.tsv'
    table.write_text('A\tX\r\n\nB\tX\nZ\tY\n')  # Windows line ending, blank line, duplicate new name, unknown name
    out = tmp_path / 'renamed.nwk'
    rename_main(['-i', str(tree_file), '-o', str(out), '-r', str(table)])
    assert tip_names(out) == ['C', 'D', 'X', 'X']
    assert 'not found in the tree: Z' in caplog.text
    assert 'shared by several tips: X' in caplog.text


def test_tree_renamer_bad_table(tree_file, tmp_path, caplog):
    table = tmp_path / 'rename.tsv'
    table.write_text('A\tX\nB X\n')
    out = tmp_path / 'renamed.nwk'
    assert exit_code(rename_main, ['-i', str(tree_file), '-o', str(out), '-r', str(table)]) == 1
    assert 'Line 2' in caplog.text
    assert not out.exists()


def test_read_rename_table_keeps_spaces(tmp_path):
    table = tmp_path / 'rename.tsv'
    table.write_text('A\tnew name \n')
    assert read_rename_table(table) == {'A': 'new name '}

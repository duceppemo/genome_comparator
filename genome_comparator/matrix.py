"""Read, validate and write square distance matrices."""

from pathlib import Path

import numpy as np
import pandas as pd


class MatrixError(Exception):
    pass


def read_matrix(path):
    """
    Read a square distance matrix (.tsv, .csv, .xlsx or .xls).
    First row and first column hold the sample names.

    :return: pandas DataFrame with identical, sorted row and column labels
    """
    path = Path(path)
    ext = path.suffix.lower()
    if ext in ('.xlsx', '.xls'):
        df = pd.read_excel(path, index_col=0, header=0)
    elif ext == '.csv':
        df = pd.read_csv(path, index_col=0, header=0)
    elif ext in ('.tsv', '.txt', '.tab'):
        df = pd.read_csv(path, index_col=0, header=0, sep='\t')
    else:
        raise MatrixError('Invalid input file type "{}". The distance matrix must be in Excel '
                          '(".xlsx" or ".xls") or text format (".csv" or ".tsv")'.format(ext))
    return validate(df)


def validate(df):
    """Check that a distance matrix is square, symmetric and complete, and sort it by sample name."""
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    if df.shape[0] != df.shape[1]:
        raise MatrixError('Distance matrix is not square ({} rows x {} columns)'.format(*df.shape))
    if df.index.has_duplicates or df.columns.has_duplicates:
        raise MatrixError('Duplicated sample names in distance matrix')
    if set(df.index) != set(df.columns):
        missing = sorted(set(df.index) ^ set(df.columns))
        raise MatrixError('Row and column names differ: {}'.format(', '.join(missing[:10])))

    labels = sorted(df.columns)
    df = df.loc[labels, labels].apply(pd.to_numeric, errors='coerce')
    values = df.to_numpy(dtype=float)

    if np.isnan(values).any():
        raise MatrixError('Distance matrix has missing or non-numeric values')
    if not np.allclose(values, values.T, atol=1e-6):
        raise MatrixError('Distance matrix is not symmetric')
    if (values < 0).any():
        raise MatrixError('Distance matrix has negative values')

    # Remove rounding noise so downstream tools that require exact symmetry are happy
    values = (values + values.T) / 2
    np.fill_diagonal(values, 0)
    return pd.DataFrame(values, index=labels, columns=labels)


def write_tsv(df, path):
    """Write a square matrix. The "#query" corner cell matches the "mash dist -t" format."""
    df.to_csv(path, sep='\t', index_label='#query', float_format='%.6g')


def write_phylip(df, path):
    """Write a square matrix in relaxed phylip format (for rapidnj, fastme, quicktree, etc.)."""
    with open(path, 'w') as f:
        f.write('{}\n'.format(len(df)))
        for name, row in zip(df.index, df.to_numpy()):
            f.write('{}\t{}\n'.format(name.replace(' ', '_'), '\t'.join('{:.6g}'.format(v) for v in row)))

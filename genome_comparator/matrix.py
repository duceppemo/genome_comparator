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
    # Everything is read as text so sample names such as "001" or "1e3" are kept as written. The first column
    # becomes the index afterwards: older pandas versions ignore dtype for the index column.
    if ext in ('.xlsx', '.xls'):
        try:
            df = pd.read_excel(path, header=0, dtype=str)
        except ImportError:
            raise MatrixError('Reading "{}" files requires an extra package: install it with '
                              '"conda install -c conda-forge xlrd", or save the file as .xlsx or .tsv'.format(ext))
    elif ext == '.csv':
        df = pd.read_csv(path, header=0, dtype=str)
    elif ext in ('.tsv', '.txt', '.tab'):
        df = pd.read_csv(path, header=0, sep='\t', dtype=str)
    else:
        raise MatrixError('Invalid input file type "{}". The distance matrix must be in Excel '
                          '(".xlsx" or ".xls") or text format (".csv" or ".tsv")'.format(ext))
    if df.shape[1] == 0:
        raise MatrixError('Distance matrix is empty')
    return validate(df.set_index(df.columns[0]))


def labels_to_str(labels):
    """Sample names as strings. Integral numbers (e.g. numeric IDs read from Excel) are written without ".0"."""
    return [str(int(x)) if isinstance(x, float) and x.is_integer() else str(x) for x in labels]


def validate(df):
    """Check that a distance matrix is square, symmetric and complete, and sort it by sample name."""
    df.index = labels_to_str(df.index)
    df.columns = labels_to_str(df.columns)

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

    if not np.isfinite(values).all():
        raise MatrixError('Distance matrix has missing, infinite or non-numeric values')
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

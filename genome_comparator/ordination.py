"""Principal coordinates analysis (PCoA) of a distance matrix with an interactive plot."""

import warnings

import pandas as pd
import plotly.express as px
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa as skbio_pcoa

# Above this many samples, use the fast approximate method instead of a full eigendecomposition
FSVD_THRESHOLD = 1000


def pcoa(df, dimensions=3):
    """
    PCoA (a.k.a. classical multidimensional scaling) is the PCA equivalent for distance matrices.

    :return: (coordinates DataFrame indexed by sample, list of % variance explained per axis)
    """
    n = len(df)
    dimensions = min(dimensions, n - 1)
    method = 'fsvd' if n > FSVD_THRESHOLD else 'eigh'
    dm = DistanceMatrix(df.to_numpy(), list(df.index))
    with warnings.catch_warnings():
        # Negative eigenvalues are expected with non-Euclidean distances such as Mash's
        warnings.simplefilter('ignore', RuntimeWarning)
        result = skbio_pcoa(dm, method=method, dimensions=dimensions)

    coords = result.samples.iloc[:, :dimensions].copy()
    coords.index = list(df.index)
    coords.columns = ['PC{}'.format(i + 1) for i in range(dimensions)]
    explained = [float(x) * 100 for x in result.proportion_explained.iloc[:dimensions]]
    return coords, explained


def read_metadata(path, color_by=None):
    """Read a tab-separated metadata file. First column must hold the sample names."""
    meta = pd.read_csv(path, sep='\t', index_col=0, dtype=str)
    meta.index = meta.index.astype(str)
    if color_by and color_by not in meta.columns:
        raise ValueError('Column "{}" not found in metadata file. Available columns: {}'.format(
            color_by, ', '.join(meta.columns)))
    return meta


def plot_pcoa(coords, explained, html_file, metadata=None, color_by=None, title=None):
    """
    Save an interactive scatter plot (self-contained html, works offline).
    Points can be coloured by a metadata column.
    """
    data = coords.copy()
    data.index.name = 'sample'
    data = data.reset_index()
    hover = list()
    if metadata is not None:
        data = data.merge(metadata, how='left', left_on='sample', right_index=True)
        hover = [c for c in metadata.columns if c != color_by]

    labels = {c: '{} ({:.1f}%)'.format(c, e) for c, e in zip(coords.columns, explained)}
    fig = px.scatter(data, x='PC1', y='PC2',
                     color=color_by, hover_name='sample', hover_data=hover,
                     labels=labels, title=title)
    fig.update_traces(marker_size=7)
    fig.write_html(html_file, include_plotlyjs=True)

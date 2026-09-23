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


# Okabe-Ito colours, in a fixed order. Checked for colour vision deficiencies on a white background with every pair
# of colours next to each other (as in a scatter plot): the first 5 are distinguishable by colour alone, the 6th only
# with a second cue, so every category also gets its own marker shape.
PALETTE = ('#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7')
SYMBOLS = ('circle', 'square', 'diamond', 'triangle-up', 'cross', 'triangle-down', 'star')
NEUTRAL = '#8C8C8C'  # "Unknown" (sample missing from the metadata) and "Other" (categories past the last style)
SINGLE_COLOUR = '#0072B2'  # All points when there is no --color-by


def category_styles(values):
    """
    Assign a unique (colour, symbol) pair to each category. Colours and symbols cycle at different lengths (6 and 7),
    giving 42 unique pairs. Categories are sorted by name so a category keeps its style between runs.
    When there are more categories, the least frequent ones are grouped into "Other".

    :param values: pandas Series of category labels (missing values = sample absent from the metadata)
    :return: (Series of category labels to plot, dict {label: (colour, symbol)}, list of labels in legend order)
    """
    values = values.fillna('Unknown').astype(str)
    categories = sorted(set(values) - {'Unknown'})
    max_styles = len(PALETTE) * len(SYMBOLS)
    if len(categories) > max_styles:
        keep = set(values[values != 'Unknown'].value_counts().index[:max_styles - 1])
        values = values.where(values.isin(keep) | (values == 'Unknown'), 'Other')
        categories = sorted(keep)
    styles = {c: (PALETTE[i % len(PALETTE)], SYMBOLS[i % len(SYMBOLS)]) for i, c in enumerate(categories)}
    order = list(categories)
    for extra in ('Other', 'Unknown'):
        if (values == extra).any():
            styles[extra] = (NEUTRAL, 'circle-open')
            order.append(extra)
    return values, styles, order


def pcoa_figure(coords, explained, metadata=None, color_by=None, title=None):
    """
    Interactive PCoA scatter plot (PC1 vs PC2). Hovering over a point shows the sample name and its metadata.
    Points can be coloured by a metadata column, with a colourblind-friendly palette and one marker shape per category.
    """
    data = coords.copy()
    data.index.name = 'sample'
    data = data.reset_index()
    hover = list()
    if metadata is not None:
        data = data.merge(metadata, how='left', left_on='sample', right_index=True)
        hover = [c for c in metadata.columns if c != color_by]

    labels = {c: '{} ({:.1f}%)'.format(c, e) for c, e in zip(coords.columns, explained)}
    style = dict()
    if color_by:
        data[color_by], styles, order = category_styles(data[color_by])
        style = dict(color=color_by, symbol=color_by, category_orders={color_by: order},
                     color_discrete_map={c: s[0] for c, s in styles.items()},
                     symbol_map={c: s[1] for c, s in styles.items()})
    else:
        style = dict(color_discrete_sequence=[SINGLE_COLOUR])

    fig = px.scatter(data, x='PC1', y='PC2', hover_name='sample', hover_data=hover, labels=labels, title=title,
                     template='plotly_white', **style)
    # White outline keeps overlapping points apart
    fig.update_traces(marker=dict(size=10, line=dict(width=1, color='white')))
    return fig


def plot_pcoa(coords, explained, html_file, metadata=None, color_by=None, title=None):
    """Save the interactive PCoA plot as a self-contained html file (works offline)."""
    pcoa_figure(coords, explained, metadata, color_by, title).write_html(html_file, include_plotlyjs=True)

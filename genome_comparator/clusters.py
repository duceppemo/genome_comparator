"""Group samples into clusters at distance thresholds (single linkage)."""

import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


def column_name(threshold):
    """"cluster_0.001" for a threshold of 0.001."""
    return 'cluster_{:g}'.format(threshold)


def single_linkage_clusters(df, threshold):
    """
    Single linkage clusters: two samples are in the same cluster if they are linked by a chain of samples, each at
    most `threshold` from the next. Unlike a cut of the UPGMA tree, the result never depends on sample order or ties.

    Clusters are numbered from 1 by decreasing size, ties broken by the first sample name, so the numbering is stable
    between runs on the same data. Samples with no neighbour within the threshold get their own cluster.

    :param df: square distance matrix (pandas DataFrame)
    :return: pandas Series {sample: cluster number}
    """
    linked = csr_matrix(df.to_numpy() <= threshold)
    _, labels = connected_components(linked, directed=False)
    members = pd.Series(df.index, index=df.index).groupby(labels).agg(list)
    order = sorted(members, key=lambda names: (-len(names), min(names)))
    number = {name: i for i, names in enumerate(order, 1) for name in names}
    return pd.Series([number[name] for name in df.index], index=df.index, name=column_name(threshold))


def cluster_table(df, thresholds):
    """
    One column of cluster numbers per threshold, from the smallest to the largest threshold.

    :return: pandas DataFrame indexed by sample name
    """
    return pd.concat([single_linkage_clusters(df, t) for t in sorted(set(thresholds))], axis=1)


def summary(column):
    """Short description of one clustering, e.g. "3 clusters, 1 with several samples (largest: 5 samples)"."""
    sizes = column.value_counts()
    return '{} cluster(s), {} with several samples (largest: {} sample(s))'.format(
        len(sizes), (sizes > 1).sum(), sizes.max())

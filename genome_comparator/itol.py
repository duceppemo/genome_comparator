"""Annotation files for the iTOL tree viewer (https://itol.embl.de): one colour strip per metadata column."""

import re

from .ordination import PALETTE, category_styles

LEGEND_SQUARE = 1  # iTOL legend shape


def colour_strip(values, title):
    """
    iTOL DATASET_COLORSTRIP file for one column. Categories get the same colours as in the PCoA plot. Past 6
    categories colours repeat, so each strip also shows its category as text.
    Samples with no value (missing from the metadata, or alone in their cluster) are left blank.

    :param values: pandas Series {sample: category}, missing values allowed
    """
    values = values.dropna().astype(str)
    labels, styles, order = category_styles(values)
    lines = ['DATASET_COLORSTRIP',
             'SEPARATOR TAB',
             'DATASET_LABEL\t{}'.format(title),
             'COLOR\t#000000',
             'SHOW_STRIP_LABELS\t{}'.format(int(len(order) > len(PALETTE))),
             'STRIP_LABEL_COLOR\t#000000',
             'LEGEND_TITLE\t{}'.format(title),
             'LEGEND_SHAPES\t' + '\t'.join([str(LEGEND_SQUARE)] * len(order)),
             'LEGEND_COLORS\t' + '\t'.join(styles[c][0] for c in order),
             'LEGEND_LABELS\t' + '\t'.join(order),
             'DATA']
    lines += ['{}\t{}\t{}'.format(sample, styles[label][0], label) for sample, label in labels.items()]
    return '\n'.join(lines) + '\n'


def worth_a_strip(values):
    """
    False for columns that would not group anything: empty, or a different value for every sample (e.g. strain
    names or IDs).
    """
    values = values.dropna()
    return values.duplicated().any()


def file_name(prefix, column):
    """Output file name, with characters other than letters, digits, ".", "-" and "_" replaced by "_"."""
    return '{}_itol_{}.txt'.format(prefix, re.sub(r'[^\w.-]', '_', column))

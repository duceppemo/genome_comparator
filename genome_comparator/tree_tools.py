"""Post-process Newick trees: collapse near-identical clades, rename tips."""

from collections import Counter


def mean_tip_distances(tree):
    """
    Average distance from each node to the tips below it, computed in a single post-order pass.

    :return: dict {node: average distance to its tips}
    """
    count = dict()
    total = dict()  # Sum of distances from the node to each of its tips
    for node in tree.postorder(include_self=True):
        if node.is_tip():
            count[node], total[node] = 1, 0.0
        else:
            count[node] = sum(count[c] for c in node.children)
            total[node] = sum(total[c] + count[c] * (c.length or 0.0) for c in node.children)
    return {node: total[node] / count[node] for node in count}


def collapse(tree, min_dist):
    """
    Collapse every clade whose average distance to its tips is smaller than min_dist.
    The collapsed clade becomes a tip named "<first tip> {<other tips>}".

    :return: number of collapsed clades
    """
    avg = mean_tip_distances(tree)
    collapsed = 0
    stack = list(tree.children)  # The root itself is never collapsed
    while stack:
        node = stack.pop()
        if node.is_tip():
            continue
        if avg[node] < min_dist:
            tips = [t.name for t in node.tips()]
            node.name = '{} {{{}}}'.format(tips[0], ','.join(tips[1:]))
            for child in list(node.children):
                node.remove(child)
            collapsed += 1
        else:
            stack.extend(node.children)
    return collapsed


def read_rename_table(path):
    """Read a two-column tab-separated file: current name, new name."""
    rename = dict()
    with open(path) as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) != 2:
                raise ValueError('Line {} of "{}" does not have exactly 2 tab-separated columns'.format(n, path))
            rename[fields[0]] = fields[1]
    return rename


def rename_tips(tree, rename):
    """
    Rename tips using exact name matches (no partial matches, so "S1" never touches "S10").

    :return: (set of names from the table that were not found in the tree,
              set of names shared by several tips after renaming)
    """
    found = set()
    for tip in tree.tips():
        if tip.name in rename:
            found.add(tip.name)
            tip.name = rename[tip.name]
    duplicates = {name for name, count in Counter(tip.name for tip in tree.tips()).items() if count > 1}
    return set(rename) - found, duplicates

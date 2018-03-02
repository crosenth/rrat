#!/usr/bin/env python3
"""
Calculates median copy number correction using the University
of Michigan rrndb database.  Median copy numbers
are calculated from the lowest rank to `root`.  Empty nodes
inherit their immediate parent's median copy number.

https://rrndb.umms.med.umich.edu/static/download/
"""
import argparse
import csv
import io
import itertools
import logging
import pkg_resources
import urllib.request
import statistics
import sys
import tarfile
import zipfile

RRNDB = 'https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.4.tsv.zip'
NCBI = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'

RANK_ORDER = [
    'forma',
    'varietas',
    'subspecies',
    'species',
    'species subgroup',
    'species group',
    'subgenus',
    'genus',
    'subtribe',
    'tribe',
    'subfamily',
    'family',
    'superfamily',
    'parvorder',
    'infraorder',
    'suborder',
    'order',
    'superorder',
    'cohort',
    'infraclass',
    'subclass',
    'class',
    'superclass',
    'subphylum',
    'phylum',
    'superphylum',
    'subkingdom',
    'kingdom',
    'superkingdom',
    'root',
]


def add_arguments(parser):
    parser.add_argument(
        'rrndb',
        nargs='?',
        metavar='tsv',
        help='copy number data with columns '
             '"NCBI tax id,16S gene count"[download from rrndb]')

    parser.add_argument(
        '-V', '--version',
        action='version',
        version=pkg_resources.get_distribution('rrat').version,
        help='Print the version number and exit.')
    log_parser = parser.add_argument_group(title='logging options')
    log_parser.add_argument(
        '-l', '--log',
        metavar='FILE',
        type=argparse.FileType('a'),
        default=sys.stdout,
        help='Send logging to a file')
    log_parser.add_argument(
        '-v', '--verbose',
        action='count',
        dest='verbosity',
        default=0,
        help='Increase verbosity of screen output '
             '(eg, -v is verbose, -vv more so)')
    log_parser.add_argument(
        '-q', '--quiet',
        action='store_const',
        dest='verbosity',
        const=0,
        help='Suppress output')

    parser.add_argument(
        '--nodes',
        type=argparse.FileType('r'),
        help='location of csv nodes file with columns'
             'tax_id,parent_id,rank [download from ncbi]')
    parser.add_argument(
        '--merged',
        type=argparse.FileType('r'),
        help='location of csv merged file with columns'
             'old_tax_id,tax_id [download from ncbi]')
    parser.add_argument(
        '--out',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='output 16s rrndb with taxids and counts')

    return parser


def fix_rows(rows):
    """
    concat row pieces that are split with newlines
    """
    for r in rows:
        while len(r) != 19:
            next_row = next(rows)
            r[-1] += next_row[0]
            r.extend(next_row[1:])
        yield r


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description=__doc__)
    parser = add_arguments(parser)
    args = parser.parse_args(args)
    setup_logging(args)

    if not (args.nodes and args.merged):
        logging.info('downloading ' + NCBI)
        tar, headers = urllib.request.urlretrieve(NCBI, 'taxdump.tar.gz')
        logging.debug(str(headers).strip())
        taxdmp = tarfile.open(name=tar, mode='r:gz')

    if args.nodes:
        nodes = csv.reader(args.nodes)
    else:
        nodes = io.TextIOWrapper(taxdmp.extractfile('nodes.dmp'))
        nodes = (n.strip().replace('\t', '').split('|') for n in nodes)

    # designate None and 'root' as our traversing stop values in the tax tree
    nodes = (n[:3] for n in nodes)
    nodes = [['1', None, 'root'] if n[0] == '1' else n for n in nodes]
    parents = {n[0]: n[1] for n in nodes}
    ranks = {n[0]: n[2] for n in nodes}

    if args.merged:
        merged = csv.reader(args.merged)
    else:
        merged = io.TextIOWrapper(taxdmp.extractfile('merged.dmp'))
        merged = (m.strip().replace('\t', '').split('|') for m in merged)
    merged = dict(m[:2] for m in merged)

    if args.rrndb:
        rrndb = (row for row in csv.reader(open(args.rrndb)))
    else:
        logging.info('downloading ' + RRNDB)
        zp, headers = urllib.request.urlretrieve(RRNDB, 'rrnDB-5.4.tsv.zip')
        logging.debug(str(headers).strip())
        rrndb = io.TextIOWrapper(zipfile.ZipFile(zp).open('rrnDB-5.4.tsv'))
        rrndb = (r.strip().split('\t') for r in rrndb)
        rrndb = fix_rows(rrndb)  # remove random newlines in rows

    header = next(rrndb)
    tax_id = header.index('NCBI tax id')
    count = header.index('16S gene count')
    rrndb = ([row[tax_id], float(row[count])] for row in rrndb if row[count])
    logging.info('updating merged tax_ids')
    rrndb = ([merged.get(t, t), c] for t, c in rrndb)

    '''
    Occasionally, circular lineages can exist in ncbi lineages.  To account
    for this any child node with same rank as parent is set to no_rank and
    expanded later.
    '''
    new_nodes = []
    for i, pi, rank in nodes:
        if rank != 'no rank' and rank != 'root' and ranks[pi] == rank:
            logging.warning(i + ' has same rank as parent')
            rank = 'no rank'
            ranks[i] = rank
        new_nodes.append((i, pi, rank))
    nodes = new_nodes

    logging.info('crawling no rank lineages')
    # create mini lineages for expanding no_rank nodes
    lineages = []
    no_ranks = [n for n in nodes if n[2] == 'no rank']
    for i, _, rank in no_ranks:
        while rank == 'no rank':
            lineages.append(i)
            i = parents[i]
            rank = ranks[i]
        lineages.append(i)

    '''
    We will utilize the available rank data to help order the nodes.
    `no_rank`  will need to be expanded and inserted into the rank_order.
    '''
    rank_order = RANK_ORDER
    logging.info('expanding no ranks')
    prev_rank = ranks[lineages.pop()]  # start with last node which is ranked
    for i in reversed(lineages):
        rank = ranks[i]
        if rank == 'no rank':
            rank = prev_rank + '_'
            if rank not in rank_order:
                rank_order.insert(rank_order.index(prev_rank), rank)
            ranks[i] = rank
        prev_rank = rank

    logging.info('calculating medians')
    copy_nums = {}

    rrndb = sorted(rrndb, key=lambda x: x[0])  # key = tax_id
    for i, grp in itertools.groupby(rrndb, lambda x: x[0]):
        copy_nums[i] = statistics.median(g[1] for g in grp)

    # Calculate copy_nums starting from the bottom (forma) working up to root
    for r in rank_order:
        tax_ids = [i for i in copy_nums if ranks[i] == r]
        tax_ids = sorted(tax_ids, key=lambda x: parents[x])
        for k, g in itertools.groupby(tax_ids, key=lambda x: parents[x]):
            if k is not None:  # last node is root with a parent of None
                copy_nums[k] = statistics.median(copy_nums[i] for i in g)

    # Starting from root, traverse down.  Empty nodes will inherit parent value
    tax_ids = [i for i, _, _ in nodes if i not in copy_nums]
    rank_order = list(reversed(rank_order))  # make root rank first
    tax_ids = sorted(tax_ids, key=lambda x: rank_order.index(ranks[x]))
    for i in tax_ids:
        copy_nums[i] = copy_nums[parents[i]]

    def by_rank(row):
        return rank_order.index(ranks[row[0]])
    rows = sorted(copy_nums.items(), key=by_rank)
    writer = csv.writer(args.out)
    writer.writerow(['tax_id', 'median'])
    writer.writerows(rows)


def setup_logging(namespace):
    loglevel = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }.get(namespace.verbosity, logging.DEBUG)
    if namespace.verbosity > 1:
        logformat = '%(levelname)s rrat %(message)s'
    else:
        logformat = 'rrat %(message)s'
    logging.basicConfig(stream=namespace.log, format=logformat, level=loglevel)


if __name__ == '__main__':
    sys.exit(main())

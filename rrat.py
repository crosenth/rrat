#!/usr/bin/env python3
"""
Calculates median copy number corrections for taxonomic nodes
using the University of Michigan rrndb database.  Median copy numbers
are calculated from the lowest rank to `root`.  Empty nodes
inherit their immediate parent's median copy number.

rrndb - https://rrndb.umms.med.umich.edu/static/download/
ncbi - ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
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


def add_arguments(parser):
    parser.add_argument(
        'rrndb',
        nargs='?',
        metavar='tsv',
        help='copy number data with columns '
             '"NCBI tax id,16S gene count" [download from rrndb]')

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
        help='location of header-less csv nodes file with columns '
             'tax_id,parent_id,rank [download from ncbi]')
    parser.add_argument(
        '--merged',
        type=argparse.FileType('r'),
        help='location of header-less csv merged file with columns '
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

    # 0 tax_id, 1 parent, 2 rank
    nodes = (n[:2] for n in nodes)

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
    rrndb = sorted(rrndb, key=lambda x: x[0])  # key = tax_id

    rrndb_medians = {}
    for i, grp in itertools.groupby(rrndb, lambda x: x[0]):  # by tax_id
        rrndb_medians[i] = statistics.median(g[1] for g in grp)

    class Node:
        def __init__(self, tax_id, rrndb):
            self.tax_id = tax_id
            self.children = []
            self.median = rrndb

        def add_child(self, child):
            self.children.append(child)

        def set_parent(self, parent):
            self.parent = parent

        def set_median(self, median):
            if self.median is None:
                self.median = median
            else:
                raise ValueError(
                    'Node {tax_id} already contains '
                    'a median copy number {median}'.format(self))

        def post_order(self):
            if self.median is None:
                medians = []
                for c in self.children:
                    med = c.post_order()
                    if med is not None:
                        medians.append(med)
                if medians:
                    self.median = statistics.median(medians)
                else:
                    self.median = None
            return self.median

        def pre_order(self):
            for c in self.children:
                if c.median is None:
                    c.set_median(self.median)
                c.pre_order()

        def __str__(self):
            return '{} --> {}'.format(self.tax_id, self.median)

        def __repr__(self):
            return '{} --> {}'.format(self.tax_id, self.median)

    node_objs = {}

    logging.info('building node tree')
    for tax_id, parent_id in nodes:
        node = Node(tax_id, rrndb_medians.get(tax_id, None))
        node_objs[tax_id] = node
        if tax_id != '1':  # root has no parent
            if parent_id in node_objs:
                parent = node_objs[parent_id]
            else:
                parent = Node(parent_id, rrndb_medians.get(tax_id, None))
                node_objs[parent_id] = parent
            parent.add_child(node)
            node.set_parent(parent)

    print(node_objs['433'].post_order())
    return

    root = node_objs['1']
    logging.info('executing postorder traversal')
    median = root.post_order()
    print(median)


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

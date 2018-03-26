#!/usr/bin/env python3
"""
Calculates median copy number corrections for taxonomic nodes
using the University of Michigan rrndb database.  Median copy numbers
are calculated using a post-order traversal to calculate median values
followed by a pre-order traversal to fill in missing nodes by inheritence
from the parent node.

rrndb - https://rrndb.umms.med.umich.edu/static/download/
ncbi - ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
"""
import argparse
import csv
import io
import logging
import os
import pkg_resources
import urllib.request
import statistics
import sys
import tarfile
import zipfile

RRNDB = 'https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.4.tsv.zip'
NCBI = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'


class Node:
    def __init__(self, tax_id, rrndb):
        self.tax_id = tax_id
        if rrndb:
            self.median = statistics.median(rrndb)
            self.how = 'R'
        else:
            self.median = None
        self.children = []

    def __repr__(self):
        '''
        prints the tax_id, median and how the number was calculated:

        R = rrndb, M = median, I = inherited
        '''
        return '{} --> {} ({})'.format(self.tax_id, self.median, self.how)

    def add_child(self, child):
        self.children.append(child)

    def post_order(self):
        '''
        traverse through the all children (subtrees) before calculating
        the median copy number value
        '''
        if self.median is None:
            child_medians = []
            for c in self.children:
                med = c.post_order()
                if med is not None:
                    child_medians.append(med)
            if child_medians:
                self.median = statistics.median(child_medians)
                self.how = 'M'
        return self.median

    def pre_order(self, print_char='|_____'):
        '''
        Fill in any missing subtree medians with the parent median value
        '''
        logging.debug(print_char + str(self))
        for c in self.children:
            if c.median is None:
                c.median = self.median
                c.how = 'I'
            c.pre_order('|  ' + print_char)

    def write_tree(self, file_obj):
        file_obj.write('{},{}\n'.format(self.tax_id, self.median))
        for c in self.children:
            c.write_tree(file_obj)


def add_arguments(parser):
    parser.add_argument(
        'copy_nums',
        nargs='?',
        metavar='zip',
        help='copy number data with columns '
             '"NCBI tax id,16S gene count"')

    parser.add_argument(
        '--root',
        default='1',
        metavar='',
        help='root node id [%(default)s]')
    parser.add_argument(
        '--nodes',
        metavar='',
        type=argparse.FileType('r'),
        help='location of header-less csv nodes file with columns '
             'tax_id,parent_id,rank [download from ncbi]')
    parser.add_argument(
        '--merged',
        metavar='',
        type=argparse.FileType('r'),
        help='location of header-less csv merged file with columns '
             'old_tax_id,tax_id [download from ncbi]')
    parser.add_argument(
        '-V', '--version',
        action='version',
        version=pkg_resources.get_distribution('rrat').version,
        help='Print the version number and exit.')

    url_parser = parser.add_argument_group(title='urls')
    url_parser.add_argument(
        '--rrndb',
        default=RRNDB,
        help='[%(default)s]')
    url_parser.add_argument(
        '--ncbi',
        default=NCBI,
        help='[%(default)s]')

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
        '--out',
        metavar='',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='output 16s rrndb with taxids and counts')

    return parser


def build_tree(nodes, medians, root):
    tree = {}
    for tax_id, parent_id in nodes:
        if tax_id in tree:
            node = tree[tax_id]
        else:
            node = Node(tax_id, medians.get(tax_id, None))
            tree[tax_id] = node

        if tax_id == root:  # root has no parent
            continue

        if parent_id in tree:
            parent = tree[parent_id]
        else:
            parent = Node(parent_id, medians.get(tax_id, None))
            tree[parent_id] = parent

        parent.add_child(node)
    return tree[root]


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
        logging.info('downloading ' + args.ncbi)
        tar, headers = urllib.request.urlretrieve(
           args.ncbi, os.path.basename(args.ncbi))
        logging.debug(str(headers).strip())
        taxdmp = tarfile.open(name=tar, mode='r:gz')

    if args.nodes:
        nodes = csv.reader(args.nodes)
    else:
        nodes = io.TextIOWrapper(taxdmp.extractfile('nodes.dmp'))
        nodes = (n.strip().replace('\t', '').split('|') for n in nodes)

    # tax_id,parent
    nodes = (n[:2] for n in nodes)

    if args.merged:
        merged = csv.reader(args.merged)
    else:
        merged = io.TextIOWrapper(taxdmp.extractfile('merged.dmp'))
        merged = (m.strip().replace('\t', '').split('|') for m in merged)
    merged = dict(m[:2] for m in merged)

    if args.copy_nums:
        zp = args.copy_nums
    else:
        logging.info('downloading ' + args.rrndb)
        zp, headers = urllib.request.urlretrieve(
            args.rrndb, os.path.basename(args.rrndb))
        logging.debug(str(headers).strip())
    tsv = os.path.basename(os.path.splitext(str(zp))[0])
    rrndb = io.TextIOWrapper(zipfile.ZipFile(zp).open(tsv))
    rrndb = (r.strip().split('\t') for r in rrndb)
    rrndb = fix_rows(rrndb)  # remove random newlines in rows
    header = next(rrndb)
    tax_id = header.index('NCBI tax id')
    count = header.index('16S gene count')
    rrndb = ([row[tax_id], float(row[count])] for row in rrndb if row[count])
    logging.info('updating merged tax_ids')
    rrndb = ([merged.get(t, t), c] for t, c in rrndb)

    # group copy numbers by tax_id
    copy_nums = {}
    for i, n in rrndb:
        if i in copy_nums:
            copy_nums[i].append(n)
        else:
            copy_nums[i] = [n]

    logging.info('building node tree')
    root = build_tree(nodes, copy_nums, args.root)
    logging.info('calculating medians by post-order traversal')
    root.post_order()
    logging.info('assigning empty nodes by pre-order traversal')
    root.pre_order()
    logging.info('writing to file')
    args.out.write('tax_id,median\n')
    root.write_tree(args.out)


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

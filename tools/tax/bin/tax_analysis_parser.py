#!/usr/bin/env python3
import sys
import argparse
import collections
import itertools
import logging

from lxml import etree
from lxml.builder import E

import gettax

__version__ = '0.81'
logger = logging.getLogger('tax_analysis_parser')

def get_lineage(tax_id, cache, conn):
    '''returns lineage fro given tax_id'''
    res = cache.get(tax_id)
    if res is None:
        if tax_id in (0, 1):
            res = [1]
        else:
            tax_info = conn.gettax(tax_id)
            parent_tax_id = tax_info.parent_tax_id if tax_info else 0
            parent_lineage = get_lineage(parent_tax_id, cache, conn)
            res = list(parent_lineage)
            res.append(tax_id)
            cache[tax_id] = res
    return res

def get_or_add_node(nodes, counter, conn, tax_id):
    '''finds or creates taxon xml node for given tax_id'''
    node = nodes.get(tax_id)
    if node is not None:
        return node

    try:
        count = counter.pop(tax_id)
    except KeyError:
        count = 0
    node = E.taxon(tax_id=str(tax_id),
                   self_count=str(count))

    tax_info = conn.gettax(tax_id)
    if tax_info:
        if tax_info.rank:
            node.attrib['rank'] = tax_info.rank
        node.attrib['name'] = tax_info.scientific_name
        parent_tax_id = tax_info.parent_tax_id
    else:
        node.attrib['name'] = 'unknown'
        parent_tax_id = 1

    nodes[tax_id] = node
    if tax_id != 1:
        parent_node = get_or_add_node(nodes, counter, conn, parent_tax_id)
        parent_node.append(node)
    return node

def calculate_total_counts(node):
    total_count = int(node.attrib['self_count'])
    for child in node:
        calculate_total_counts(child)
        total_count += int(child.attrib['total_count'])
    node[:] = sorted(node, key=lambda child: int(child.attrib['total_count']), reverse=True)
    node.attrib['total_count'] = str(total_count)

def build_tree(counter, conn):
    '''builds xml tree'''
    logger.info('Building xml tree')
    nodes = {}
    while counter:
        for tax_id in counter: # get random key
            break
        get_or_add_node(nodes, counter, conn, tax_id)
    assert not counter, 'must be consumed'
    root = nodes[1]

    logger.info('Calculating totals')
    calculate_total_counts(root)
    return root

def deduce_tax_id(hits, lineage_cache, conn):
    '''deduces single tax_id out of set of hits'''
    assert hits
    if len(hits) == 1:
        for hit in hits: # get first (and only) element of set
            return hit

    # choose closest common ancestor
    # then move to the most specific
    # non-ambiguous tax
    #
    # example (hits are in parens):
    # (a) - b - c - (d)
    #         +- (e)
    # closest common ancestor for all hits: (a)
    # most specific without ambiguity: b
    # result is: b

    lineages = [get_lineage(hit, lineage_cache, conn) for hit in hits]
    last_matching_tax_id = None
    for nodes in itertools.zip_longest(*lineages):
        first_not_empty_node = None
        for node in nodes:
            if node is not None:
                if first_not_empty_node is None:
                    first_not_empty_node = node
                elif node != first_not_empty_node:
                    #print 'deduced', hits, last_matching_tax_id
                    return last_matching_tax_id
        last_matching_tax_id = first_not_empty_node
    #print 'deduced', hits, last_matching_tax_id
    return last_matching_tax_id

def iterate_merged_spots(f, collated):
    last_spot = None
    last_hits = None
    for line in f:
        parts = line.split('\t')
        spot = parts[0]
        hits = set(int(p.split('x')[0]) for p in parts[1:])
        assert hits, 'no hits for spot: %s' % spot
        if collated:
            yield hits
        else:
            if spot == last_spot:
                last_hits |= hits
            else:
                if last_spot:
                    assert last_spot <= spot, 'input is not sorted'
                    yield last_hits
                last_spot = spot
                last_hits = hits
    if last_spot:
        yield last_hits

def iterate_merged_spots_compact(f, collated):
    last_spot = None
    last_hits = None
    for line in f:
        if not line:
            break

        if line[0] != '\t':
            if last_spot:
                last_spot = None
                yield last_hits, 1

            parts = line.split('\t')
            copies = int(parts[0])
            hits = set(int(p.split('x')[0]) for p in parts[1:])
            yield hits, copies
        else:
            line = line[1:]
            parts = line.split('\t')
            spot = parts[0]
            hits = set(int(p.split('x')[0]) for p in parts[1:])
            assert hits, 'no hits for spot: %s' % spot
            if collated:
                yield hits, 1
            else:    
                if spot == last_spot:
                    last_hits |= hits
                else:
                    if last_spot:
                        assert last_spot <= spot, 'input is not sorted'
                        yield last_hits, 1
                    last_spot = spot
                    last_hits = hits

    if last_spot:
        yield last_hits, 1

def parse(f, conn, wgs_mode, compact, collated, include_tax_ids=[]):
    '''parses tax_analysis output file'''
    counter = collections.Counter()
    counter[1] = 0 # explicitly add root
    for tax_id in include_tax_ids:
        counter[tax_id] = 0

    lineage_cache = {}

    if wgs_mode:
        for line in f:
            parts = line.split('\t')
            hits = parts[1:]
            assert hits, 'no hits for spot: %s' % parts[0]
            for hit in hits:
                if 'x' in hit:
                    tax_id, count = map(int, hit.split('x'))
                else:
                    tax_id = int(hit)
                    count = 1
                counter[tax_id] += count
    else:
        if compact:
            for hits, copies in iterate_merged_spots_compact(f, collated):
                tax_id = deduce_tax_id(hits, lineage_cache, conn)
                counter[tax_id] += copies
        else:
            for hits in iterate_merged_spots(f, collated):
                tax_id = deduce_tax_id(hits, lineage_cache, conn)
                counter[tax_id] += 1

    xml = build_tree(counter, conn)
    return xml

def main():
    parser = argparse.ArgumentParser('tax analysis output parser, builds tree from raw output')
    parser.add_argument('-v', '--verbose', action='store_true', help='enable verbose mode')
    parser.add_argument('-t', '--tax-dump')
    parser.add_argument('-c', '--sqlite-cache')
    parser.add_argument('-r', '--rebuild-timeout', type=float, help='minimum delay between cache rebuilds in seconds')
    parser.add_argument('-i', '--include-tax-id', type=int, action='append', help='include taxon into tax tree even if it has no hits')
    parser.add_argument('--connection-timeout', type=float, help='cache connection timeout, to wait until rebuild is completed')
    parser.add_argument('--compact', action='store_true')
    parser.add_argument('--wgs-mode', action='store_true', help='''
In regular mode parser assigns single consensus tax_id for each input sequence.
It then builds hierarchy showing counts of sequences matching each particular tax_id.
It doesn't work very well for wgs, because there are only few sequences and because they can have very uneven lengths.
With this flag parser builds hierarchy based on count of kmer hits, not the count of sequences.
'''.strip())
    parser.add_argument('path', nargs='?', help='path ot file with tax analysis output, if empty reads from stdin')
    parser.add_argument('--collated', action='store_true', help='the input is collated')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    if not args.tax_dump and not args.sqlite_cache:
        parser.error('need tax_dump or sqlite_cache or both to work')

    if args.path:
        logger.info('Reading %s', args.path)
        f = open(args.path)
    else:
        logger.info('Reading stdin')
        f = sys.stdin

    with gettax.connect(args.tax_dump, args.sqlite_cache, args.rebuild_timeout, args.connection_timeout) as conn:
        xml = parse(f, conn, args.wgs_mode, args.compact, args.collated, args.include_tax_id or [])
    xml = E.taxon_tree(xml, parser_version=__version__)
    print (etree.tostring(xml, pretty_print=True).decode('utf-8'))

if __name__ == '__main__':
    main()

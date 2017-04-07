#!/opt/python-2.7env/bin/python

import sys
import argparse
import collections
import itertools
import logging
import types

from lxml import etree

logger = logging.getLogger('tax_analysis_report')

def flatten(arg):
    if isinstance(arg, list) or isinstance(arg, types.GeneratorType):
        return sum((flatten(item) for item in arg), [])
    else:
        return [arg]

def format_node(node, grand_total, offset, args):
    rank = node.attrib.get('rank')
    if not rank:
        for subnode in node:
            yield format_node(subnode, grand_total, offset, args)
    else:
        name = node.attrib['name']
        total_count = int(node.attrib['total_count'])
        if total_count < args.cutoff_hit_count:
            return
        
        rate = float(total_count) / grand_total
        percent = rate*100
        if percent < args.cutoff_percent:
            return

        if args.no_padding:
            percent_precision = ''
            hits_precision = ''
        else:
            percent_precision = '5'
            hits_precision = str(len(str(grand_total)))
            
        pattern = '%s%s\t%' + percent_precision + '.2f%%  (%' + hits_precision + 'd hits)'
        yield pattern % (offset, name, percent, total_count)
        
        for subnode in node:
            yield format_node(subnode, grand_total, offset+args.indent, args)

def pad_tree(lines, separator):
    max_name_len = 0
    for line in lines:
        name, stats = line.rsplit('\t', 1)
        max_name_len = max(max_name_len, len(name))

    for line in lines:
        name, stats = line.rsplit('\t', 1)
        yield name.ljust(max_name_len) + separator + stats

def format_tax_tree(tree, args):
    if len(tree) == 0:
        return 'tree is empty'
    root = tree[0]
    total_count = int(root.attrib['total_count'])

    res = [format_node(node, total_count, '', args) for node in root]
    res = flatten(res)
    if args.no_padding:
        for idx in range(len(res)):
            name, stats = res[idx].rsplit('\t', 1)
            res[idx] = name + args.separator + stats
    else:
        res = pad_tree(res, args.separator)
    return '\n'.join(res)

def main():
    parser = argparse.ArgumentParser(description='tax analysis report creation tool, builds human-readable report from xml tree')
    parser.add_argument('-v', '--verbose', action='store_true', help='enable verbose mode')
    formatting = parser.add_argument_group('formatting options')
    formatting.add_argument('--indent', metavar='STR', help='indentation string, default is two spaces', default='  ')
    formatting.add_argument('--separator', metavar='STR', help='name/stats separator string, default is four spaces', default='    ')
    formatting.add_argument('--no-padding', action='store_true', help='disable tree padding')
    filtering = parser.add_argument_group('filtering options')
    filtering.add_argument('--cutoff-percent', metavar='NUMBER', help='cutoff percent, default is %(default)s', default='0.01', type=float)
    filtering.add_argument('--cutoff-hit-count', metavar='NUMBER', help='cutoff hit count, disabled by default', default='0', type=int)
    parser.add_argument('path', nargs='?', help='path ot file with tax analysis xml tree, if not specified reads from stdin')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)
    
    if args.path:
        logger.info('Reading %s', args.path)
        xml = etree.parse(args.path)
    else:
        logger.info('Reading stdin')
        xml = etree.parse(sys.stdin)

    print format_tax_tree(xml.getroot(), args)

if __name__ == '__main__':
    main()

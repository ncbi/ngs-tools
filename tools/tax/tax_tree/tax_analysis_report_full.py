#!/opt/python-2.7/bin/python
import sys
import collections
import argparse
import taxonomy
import report_tax_tree_full

def load_dbs_out(f):
	tax_counts = collections.Counter()

	while True:
		seq_id = f.readline()
		if not seq_id:
			return tax_counts
		
		hits = f.readline()
		assert hits, 'Odd count of lines in input'
		
		for hit in hits.split():
			tax_id = int(hit.split('x')[0])
			tax_counts[tax_id] += 1

def main():
	parser = argparse.ArgumentParser(description='tax analysis report generator')
	parser.add_argument('-l', '--lineage', action='store_true', help='build lineage')
	parser.add_argument('-t', '--tree', action='store_true', help='build tree')
	parser.add_argument('--min-count', metavar='N', type=int, help='cut items with less than given amount of hits')
	parser.add_argument('--min-items', metavar='N', type=int, help='minimum amount of items to show before cutting, default is %(default)s', default=10)
	parser.add_argument('path', help='path to dbs output, pass "-" to read from stdin')

	args = parser.parse_args()
	if args.path == '-':
		f = sys.stdin
	else:
		f = open(args.path)
	tax_counts = load_dbs_out(f)
	f.close()

	tax_counts = sorted(tax_counts.items(), key=lambda x: x[1], reverse=True)
	lineage_count = []
	for idx, tax_pair in enumerate(tax_counts):
		tax_id, count = tax_pair
		if args.min_count is not None and count < args.min_count and idx > args.min_items:
			print "..."
			break

		output = '%s x %s' % (count, tax_id)
		if args.lineage or args.tree:
			lineage = taxonomy.get_tax_lineage(tax_id)
			lineage_count.append((lineage, count))
			if args.lineage:
				output += ' ' + str(lineage)
		print output

	if args.tree:
		print ""
		print ""
		
		tree, weight = report_tax_tree_full.build_tree(lineage_count)
		report_tax_tree_full.print_tree(tree)
	
if __name__ == "__main__":
	main()

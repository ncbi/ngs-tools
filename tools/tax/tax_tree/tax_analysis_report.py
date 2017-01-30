#!/opt/python-2.7/bin/python
import sys
import collections
import argparse
import taxonomy
import report_tax_tree

def load_dbs_out(f):
	tax_counts = collections.Counter()
	read_and_taxes = []

	while True:
#		seq_id = f.readline()
#		if not seq_id:
#			return tax_counts, read_and_taxes
		
		hits = f.readline()
		if not hits:
			return tax_counts, read_and_taxes

		hits = hits.split()

		seq_id = hits[0]
		hits = hits[1:]

		assert hits, 'ho hits for seq_id' + seq_id
		
		taxes = []
		for hit in hits:
			tax_id = int(hit.split('x')[0])
			tax_counts[tax_id] += 1
			taxes.append(tax_id)

		read_and_taxes.append((seq_id, taxes))

def print_tax_counts(tax_counts, args):
	tax_counts = sorted(tax_counts.items(), key=lambda x: x[1], reverse = True)
	for idx, tax_pair in enumerate(tax_counts):
		tax_id, count = tax_pair
		if args.min_count is not None and count < args.min_count and idx > args.min_items:
			print "..."
			break

		output = '%s x %s' % (report_tax_tree.to_str(count), tax_id)
		if args.lineage:
			lineage = get_tax_lineage(tax_id)
			output += ' ' + str(lineage)
		print output

def get_lineage_count(tax_counts, args):
	tax_counts = sorted(tax_counts.items(), key=lambda x: x[1], reverse = True)
	lineage_count = []
	for idx, tax_pair in enumerate(tax_counts):
		tax_id, count = tax_pair
#		if args.min_count is not None and count < args.min_count and idx > args.min_items:
#			break

		lineage_count.append((get_tax_lineage(tax_id), count))

	return lineage_count

lineages_by_tax_id = dict()

def get_tax_lineage(tax_id):
	global lineages_by_tax_id

	if not tax_id in lineages_by_tax_id:
		lineages_by_tax_id[tax_id] = taxonomy.get_tax_lineage(tax_id)

	return lineages_by_tax_id[tax_id]

def consensus_tax(current_taxes, tax_hits):
	pop_tax_id = None
	pop_tax_id_count = 0

	for tax_id in current_taxes:
		tax_count = tax_hits[tax_id]
		if tax_count > pop_tax_id_count:
			pop_tax_id = tax_id
			pop_tax_id_count = tax_count

	return pop_tax_id

def spot_of_read(read):
	dot_pos = read.find('.')
	if dot_pos < 0:
		return read

	return int(read[:dot_pos]) # todo: remove int ?

def get_tax_spots(tax_hits, read_and_taxes):
	tax_counts = collections.Counter()

	current_spot = None
	current_taxes = set()
	for read, taxes in read_and_taxes:
		spot = spot_of_read(read)
		if current_spot != spot:
			if current_spot != None:
				tax_counts[consensus_tax(current_taxes, tax_hits)] += 1
			current_spot = spot
			current_taxes = set()

		current_taxes.update(taxes)

	if current_taxes:
		tax_counts[consensus_tax(current_taxes, tax_hits)] += 1

	return tax_counts

def sum_count(tax_counts):
	s = 0
	for tax in tax_counts:
		s += tax_counts[tax]

	return s

def convert_to_percent(tax_counts, total):
	if not total:
		return tax_counts

	total = total * 1.0
	for tax in tax_counts:
		tax_counts[tax] = 100.0 * tax_counts[tax]/total

	return tax_counts

def main():
	parser = argparse.ArgumentParser(description='tax analysis report generator')
	parser.add_argument('-l', '--lineage', action='store_true', help='build lineage')
	parser.add_argument('-t', '--tree', action='store_true', help='build tree')
	parser.add_argument('-p', '--percent', action='store_true', help='percent instead of count', default = False)
#	parser.add_argument('-hs', '--hits', action='store_true', help='print taxonomy hit statistics')
	parser.add_argument('--min-count', metavar='N', type=int, help='cut items with less than given amount of hits', default = 4)
	parser.add_argument('--min-items', metavar='N', type=int, help='minimum amount of items to show before cutting, default is %(default)s', default = 10)
	parser.add_argument('path', help='path to dbs output, pass "-" to read from stdin')

	args = parser.parse_args()
	if args.path == '-':
		f = sys.stdin
	else:
		f = open(args.path)

	tax_hits, read_and_taxes = load_dbs_out(f)
	f.close()

	tax_spots = get_tax_spots(tax_hits, read_and_taxes)
	if args.percent:
		s = sum_count(tax_spots)
#		print "sum count is", s
		tax_spots = convert_to_percent(tax_spots, s)
		
	print_tax_counts(tax_spots, args)

#	if args.hits: # old style statistics
#		print ""
#		print ""
#		print_tax_counts(tax_hits, args)

	if args.tree:
		print ""
		print ""
		
		tree, weight = report_tax_tree.build_tree(get_lineage_count(tax_spots, args))
		report_tax_tree.print_tree(tree, 0, args.min_count)
	
if __name__ == "__main__":
	main()

#!/usr/bin/python
import sys
import taxonomy
import report_tax_tree

def taxid_count(line):
	line = line.split('x')
	if len(line) == 1:
		return int(line[0]), 1
	else:
		return int(line[0]), int(line[1])

def parse_taxes(line):
	line = line.split()
	
	if not line:
		print line
		raise "something wrong with the line"

	return [taxid_count(l) for l in line]

lineage_by_tax = dict()
def get_tax_lineage(tax_id):
	if not tax_id in lineage_by_tax:
		lineage_by_tax[tax_id] = taxonomy.get_tax_lineage(tax_id)

	return lineage_by_tax[tax_id]

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <dbs out>"
		return

	f = open(sys.argv[1])
	while True:
		seq_id = f.readline().rstrip()
		if not seq_id:
			return
		taxes = parse_taxes(f.readline())
		print seq_id
		for tax, count in taxes:
			print count, 'x', get_tax_lineage(tax)

		print ""

main()

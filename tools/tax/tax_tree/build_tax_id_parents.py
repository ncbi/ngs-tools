#!/usr/bin/python
import sys
import taxonomy

def load_tax_ids(filename):
	taxes = set()
	f = open(filename)
	for line in f:
		line = line.split('/')[2:]
		tax_id = int(line[-1].split('.')[0])
		taxes.add(tax_id)

	return taxes

ROOT = 1
parent_of = dict()
def print_parents(tax_id):
	global parent_of

	if tax_id in parent_of:
		return

	parent = taxonomy.get_parent_tax_id_of_species(str(tax_id))
	parent_of[tax_id] = parent
	print tax_id, parent
	if parent != ROOT:
		print_parents(parent)

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <files.list>"
		return

	tax_ids = load_tax_ids(sys.argv[1])
	for tax_id in tax_ids:
		print_parents(tax_id)

main()

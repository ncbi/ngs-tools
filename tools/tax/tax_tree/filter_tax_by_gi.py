#!/usr/bin/python
import sys

def load_taxes(filename):
	taxes = set()
	f = open(filename)
	for line in f:
		line = line.split('/')[-1]
		tax = line.split(".fasta")[0]
		taxes.add(int(tax))

	return taxes

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 3:
		print "need <tax by gi file> <files.list>"
		return

	taxes_to_keep = load_taxes(sys.argv[2])
	f = open(sys.argv[1])

	for line in f:
		line = line.split()
		gi = int(line[0])
		tax_id = int(line[1])

		if tax_id in taxes_to_keep:
			print gi, tax_id
	
main()

#!/opt/python-2.7/bin/python
import sys
import argparse

def load_report_out(filename):
	f = open(filename)
	taxes = set()

	while True:
		line = f.readline()
		if not line:
			break

		hits = line.split('\t')[1:]
		for h in hits:
			tax = int(h.split('x')[0])
			taxes.add(tax)

	return taxes

def main():
	if len(sys.argv) < 2:
		print >> sys.stderr, "need <.hits file>"
		return

	taxes = load_report_out(sys.argv[1])
	for tax in taxes:
		print tax
	
if __name__ == "__main__":
	main()

#!/opt/python-2.7/bin/python
import sys
import argparse

def load_report_out(f):
	tax_counts = []

	while True:
		line = f.readline()
		if not line:
			break

		line = line.split(" x ")
		if len(line) < 2:
			break

		count = int(line[0])

		tax = int(line[1].split()[0])		

		tax_counts.append((tax, count))

	return tax_counts

def main():
	parser = argparse.ArgumentParser(description='tax analysis report tax list generator')
	parser.add_argument('--min-count', metavar='N', type=int, help='cut items with less than given amount of hits', default = 1)
	parser.add_argument('path', help='path to report, pass "-" to read from stdin')

	args = parser.parse_args()
	if args.path == '-':
		f = sys.stdin
	else:
		f = open(args.path)

	taxes = load_report_out(f)
	for tax, count in taxes:
		if count >= args.min_count:
			print tax
	
if __name__ == "__main__":
	main()

#!/opt/python-2.7/bin/python
#!/usr/bin/python
import os
import sys

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print >> sys.stderr, "need <fastacmd -T output>"
		return

	f = open(sys.argv[1])

	detected_gis = set()
	detected_taxes = set()

	count = 0
	gi = None
	while True:
		line_seq_id = f.readline()
		if not line_seq_id:
			break

		if not line_seq_id.rstrip():
			continue

		line_tax = f.readline()

		line = f.readline()
		while line.rstrip():
			line = f.readline()

#		print "line_seq_id", line_seq_id
		gi = int( line_seq_id.split(" id: gi|")[1].split('|')[0] )
		if gi <= 0:
			print >> sys.stderr, "bad gi of: ", line_seq_id

		count += 1
		if count % 1000000 == 0:
			print >> sys.stderr, '.',
			
		if gi in detected_gis:
#			print >> sys.stderr, "gi already known:", gi
			continue

		line = line_tax.split()
		if line[0] == "NCBI" and line[1] == "taxonomy" and line[2] == "id:":
			tax_id = int(line[3])
			print gi, tax_id
			detected_gis.add(gi)
			detected_taxes.add(tax_id)
		else:
			print >> sys.stderr, "unexpected taxonomy line:", line_tax

	print >> sys.stderr, "gis: ", len(detected_gis)
	print >> sys.stderr, "taxonomies: ", len(detected_taxes)

main()

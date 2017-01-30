#!/usr/bin/python
import sys

def predicted(kmer):
	pred = 0
	for i in xrange(4, len(kmer)):
		if kmer[i] == kmer[i - 1] and kmer[i - 1] == kmer[i - 2] and kmer[i - 2] == kmer[i - 3] and kmer[i - 3] == kmer[i - 4]:
			pred += 1
		elif i >= 8 and kmer[i] == kmer[i - 2] and kmer[i - 2] == kmer[i - 4] and kmer[i - 4] == kmer[i - 6] and kmer[i - 6] == kmer[i - 8]:
			pred += 1
	
	return pred

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <db file>"
		return

	f = open(sys.argv[1])
#	min_score = int(sys.argv[2])
	MIN_SCORE = 13
	CELLULAR_ORGANISMS = 131567

	count = 0
	for line in f:
#		line = line.rstrip()
		line = line.split()
		kmer, tax_id = line[0], int(line[1])
		score = predicted(kmer)
		if score >= MIN_SCORE or tax_id == CELLULAR_ORGANISMS:
			print kmer #, predicted(line)

		count += 1
		if count % 10000 == 0:
			sys.stderr.write('.')
#		else:
#			print line
	
main()

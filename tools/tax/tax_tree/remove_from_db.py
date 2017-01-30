#!/usr/bin/python
import sys

def load_file(filename):
	kmers = set()
	f = open(filename)
	for line in f:
		kmers.add(line.rstrip())

	return kmers	

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 3:
		print "need <db file> <kmer file>"
		return

	f = open(sys.argv[1])
	kmers_to_remove = load_file(sys.argv[2])

	for line in f:
		line = line.rstrip()
		kmer = line.split()[0]
		if kmer in kmers_to_remove:
			sys.stderr.write(line + '\n')
		else:
			print line
	
main()

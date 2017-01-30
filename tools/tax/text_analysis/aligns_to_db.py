#!/usr/bin/python
import sys

def any_of(db):
	for kmer in db:
		return kmer

	return None

def load_db(filename):
	f = open(filename)
	db = set()
	for line in f:
		line = line.rstrip()
		db.add(line)

	return db, len(any_of(db))

def complement_letter(ch): # // todo: optimize ?
	if ch == 'A':
		return 'T'

	if ch == 'C':
		return 'G'

	if ch == 'T':
		return 'A'

	if ch == 'G':
		return 'C'

	return ch # todo: think
	#raise "no complement letter for " + ch

def complement(seq):
#	print "direct", seq
	result = []
	for ch in seq:
		result.append(complement_letter(ch))

	rev_compl = ''.join(result)
#	print "rev_com", rev_compl

	return rev_compl

def reverse_complement(seq):
	return complement(reverse(seq))

def reverse(seq):
	return seq[::-1]

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 3:
		print "need <db file> <clean fasta file>"
		return

	f = open(sys.argv[2])
	db, kmer_len = load_db(sys.argv[1])

	line_number = 0
	desc = None
	for line in f:
		line = line.rstrip()
		if line[0] == '>':
			desc = line
			continue
#		print line_number,
		for i in xrange(len(line) - kmer_len + 1):
			kmer = line[i : i + kmer_len]
			if kmer in db: # or reverse_complement(kmer) in db: # or reverse(kmer) in db or complement(kmer) in db :
#				if desc != None:
#					print desc,
				print kmer
#				break
		
		if line_number % 256 == 0:
			sys.stderr.write('.')
		line_number += 1
	
main()

#!/usr/bin/python
import sys
import math

def load_text(filename):
	seqs = []
	f = open(filename)
	for line in f:
		seqs.append(line.rstrip())
		
	return seqs

def add_kmers(kmers, seq, kmer_len):
	for i in xrange(len(seq) - kmer_len + 1):
		kmer = seq[i : i + kmer_len]
		kmers[kmer] = kmers.get(kmer, 0) + 1

def build_kmers(seqs, kmer_len):
	kmers = dict()
	for seq in seqs:
		add_kmers(kmers, seq, kmer_len)
#		print ".",

#	print ""
	return kmers

def kmers_sum_weight(kmers):
	return sum([count for kmer, count in kmers])

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 4:
		print "need <text> <kmer len> <min coverage>"
		return

	filename = sys.argv[1]
#	if filename.endswith(".fasta"):
	seqs = load_text(filename)
#	else:
#		seqs = load_fasta(filename)
#	print "loaded", len(seqs)
	kmers = build_kmers(seqs, int(sys.argv[2]))

	MIN_COVERAGE = int(sys.argv[3])

	kmers = [(kmer, kmers[kmer]) for kmer in kmers if kmers[kmer] >= MIN_COVERAGE]
	kmers.sort(key = lambda (kmer, count) : -count)

#	sum_weight = kmers_sum_weight(kmers)
#	print sum_weight

	for kmer, count in kmers:
#		part = 1.0*count/sum_weight
#		bits = math.log(1.0/part, 2)
#		print kmer + '|' + str(count) + '|' + str( round(bits, 4) )
		print kmer + '|' + str(count)

main()

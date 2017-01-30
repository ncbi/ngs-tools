#!/usr/bin/python
import sys
import math

def load_fasta(filename):
	seqs = []
	seq = ""
	desc = ""
	f = open(filename)
	for line in f:
		if line[0] == '>':
			if seq:
				seqs.append((desc, seq))
			seq = ""
			desc = line.rstrip()
		else:
			seq += line.rstrip().upper()
		
	seqs.append((desc, seq))

	return seqs

def next_bad_char(seq, _from):
	while True:
		if _from >= len(seq):
			return len(seq)

		ch = seq[_from]
		if ch != 'A' and ch != 'C' and ch != 'T' and ch != 'G':
			return _from

		_from += 1

def seq_split(seq):
	_from = 0
	seqs = []
	MIN_SEQ_LEN = 32
	while True:
		_to = next_bad_char(seq, _from)
		valid_seq = seq[_from: _to]
		if len(valid_seq) >= MIN_SEQ_LEN:
			seqs.append(valid_seq)

		if _to == len(seq):
			break

		_from = _to + 1

	return seqs

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <fasta>"
		return

	seqs = load_fasta(sys.argv[1])
	for desc, seq in seqs:
		cleaned_seqs = seq_split(seq)
		for i in range(len(cleaned_seqs)):
#			print desc + "_" + str(i)
			print cleaned_seqs[i]

main()

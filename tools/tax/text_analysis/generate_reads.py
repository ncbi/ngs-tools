#!/usr/bin/python
import sys
import random


def load_text(filename):
	seqs = []
	f = open(filename)
	for line in f:
		seqs.append(line.rstrip())
		
	return seqs

def random_pos(seqs_len):
	return random.randint(0, seqs_len - 1 )

def get_seq(seqs, pos):
	sum_len = 0
	for seq in seqs:
		if sum_len + len(seq) > pos:
			return seq, pos - sum_len

		sum_len += len(seq)

	raise "get_seq failed"

def breakpoint_pos(seq, local_pos, real_len):
	if local_pos < 0:
		raise "local_pos < 0"

	return local_pos + real_len > len(seq)

def make_errors(s, err_per_nucl):
	for i in range(len(s)):
		if random.randint(0, err_per_nucl - 1) == 0:
			s = s[:i] + s[len(s) - i - 1] + s[i+1 : ]

	return str(s)

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 5:
		print "need <fasta file> <read len> <count> <1 error per X nucl>"
		return

	filename = sys.argv[1]
	read_len = int(sys.argv[2])
	count = int(sys.argv[3])
	err_per_nucl = int(sys.argv[4])

	seqs = load_text(filename)
	sum_seq_len = sum([len(seq) for seq in seqs])
	results_generated = 0
	while results_generated < count:
		s_pos = random_pos(sum_seq_len - read_len)
		seq, local_pos = get_seq(seqs, s_pos)

		if breakpoint_pos(seq, local_pos, read_len):
			continue 

		s = seq[local_pos : local_pos + read_len]
		s = make_errors(s, err_per_nucl)
		print s
		results_generated += 1 

main()

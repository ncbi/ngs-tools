#!/usr/bin/python
import string
import sys

translation = string.maketrans("ACTG", "TGAC")

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

def complement_(seq):
#	print "direct", seq
	result = []
	for ch in seq:
		result.append(complement_letter(ch))

	rev_compl = ''.join(result)
#	print "rev_com", rev_compl

	return rev_compl

def complement(seq):
	return string.translate(seq, translation)

def reverse_complement(seq):
	return complement(reverse(seq))

def reverse(seq):
	return seq[::-1]

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <seq>"
		return

	seq = sys.argv[1]
	print "seq      ", seq
	print "compl    ", complement(seq)
	print "rev      ", reverse(seq)
	print "rev compl", reverse_complement(seq)

main()
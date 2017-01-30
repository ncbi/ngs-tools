#!/usr/bin/python
import sys
import math

def for_every_sequence_do(filename, function):
	f = open(filename)
	for line in f:
		res = function(line.rstrip())
		for r in res:
			print r

def translate(seq):
	return [translate_frame(seq)] +	[translate_frame(seq[1:])] + [translate_frame(seq[2:])]

aminoacid = dict()

def translate_frame(seq):
	result = []
	pos = 0
	while pos + 2 < len(seq):
#		sys.stdout.write(aminoacid[seq[pos : pos + 3]])
		result.append(aminoacid[seq[pos : pos + 3]])
		pos += 3

#	if len(seq) >= 3:
#		sys.stdout.write('\n')
	return ''.join(result)

def	fill_animoacid():
	# standard code from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
	base1 =  "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
	base2 =  "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
	base3 =  "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
	result = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

	for i in range(len(result)):
		aminoacid[base1[i] + base2[i] + base3[i]] = result[i]

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <sequence list file>"
		return

	fill_animoacid()
	for_every_sequence_do(sys.argv[1], translate)

main()

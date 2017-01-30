#!/usr/bin/python
import sys

def bad(kmer):
	letters = dict()
	for ch in kmer:
		letters[ch] = letters.get(ch, 0) + 1

	letters = [letters[ch] for ch in letters]
	letters.sort(key = lambda (number) : -number)
#	if len(letters) < 3:
#		return True

	return letters[0] >= len(kmer)*7/10
	
def predicted(kmer):
	pred = 0
	step1 = 1
	step2 = 1
	for i in xrange(4, len(kmer)):
		if kmer[i] == kmer[i - 1] and kmer[i - 1] == kmer[i - 2] and kmer[i - 2] == kmer[i - 3] and kmer[i - 3] == kmer[i - 4]:
			pred += step1
			step1 = 1
		elif i >= 8 and kmer[i] == kmer[i - 2] and kmer[i - 2] == kmer[i - 4] and kmer[i - 4] == kmer[i - 6] and kmer[i - 6] == kmer[i - 8]:
			pred += step2
			step2 = 1
		else:
			step1 = 1
			step2 = 1
	
	return pred

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 3:
		print "need <db file> <min_score>"
		return

	f = open(sys.argv[1])
	min_score = int(sys.argv[2])
	count = 0
	for line in f:
		line = line.rstrip()
		score = predicted(line)
		if score >= min_score:
			print line #, predicted(line)
#		if bad(line):
		count += 1
		if count % 10000 == 0:
			sys.stderr.write('.')
#		else:
#			print line
	
main()

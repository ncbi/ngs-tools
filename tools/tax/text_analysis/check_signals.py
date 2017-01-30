#!/usr/bin/python
import sys
import seq_translate

def load_text(filename):
	seqs = []
	f = open(filename)
	for line in f:
		seqs.append(line.rstrip())
		
	return seqs

def load_signals(filename):
	f = open(filename)
	signals = dict()
	for line in f:
		line = line.rstrip().split()
		key = line[0]
		confirmations = line[1:]
		signals[key] = set(confirmations)

	return signals

def print_keys(keys):
	if not keys:
		return

	print "keys",
	for key, pos in keys:
		print key,
	print ""

def print_confirmations(confirmations):
	if not confirmations:
		return

	print "confirmations"
	for key, s in confirmations:
		print '[', key, s, ']'

	print ""

KEY_LEN = 4

def check_frame(seq, signals):
	keys = set()
	confirmations = set()

	pos = 0
	while pos < len(seq) - KEY_LEN:
		s = seq[pos : pos + KEY_LEN]

		confirmation_found = False
		for key, key_pos in keys:
			if key_pos <= pos - KEY_LEN and s in signals[key]:
				confirmations.add((key, s))
				confirmation_found = True

		if s in signals:
			keys.add((s, pos))

		if confirmation_found:
			pos += KEY_LEN
		else:
			pos += 1

	if not keys and not confirmations:
		return

	print seq
	print_keys(keys)
	print_confirmations(confirmations)

def check_signals(seq, signals):	
	frames = seq_translate.translate(seq)
#	frames = frames + seq_translate.translate(reversed(seq)) # todo: uncomment
	for frame in frames:
		check_frame(frame, signals)

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 3:
		print "need <fasta file> <signals file>"
		return

	seqs = load_text(sys.argv[1])
	signals = load_signals(sys.argv[2])
	seq_translate.fill_animoacid() # todo: constructor

	for seq in seqs:
		check_signals(seq, signals)

main()

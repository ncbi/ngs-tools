#!/usr/bin/python
import sys
import math

MIN_KEY_WORD_LEN = 12

def add_words(next, words):
	prev_word = None
	for word in words:
		if prev_word != None and len(prev_word) >= MIN_KEY_WORD_LEN:
			if not prev_word in next:
				next[prev_word] = set([word]) # todo: better
			else:
				next[prev_word].add(word)
		
		prev_word = word

def load_next_words(filename):
	next = dict()
	f = open(filename)
	for line in f:
		line = line.rstrip().split()
		add_words(next, line)
		
	return next

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <splitted text file>"
		return

	filename = sys.argv[1]

	next = load_next_words(filename)
	result = [(word, next[word]) for word in next]
	result.sort( key = lambda (word, words) : -len(words) )

	for word, words in result:
		print word,
		for w in words:
			print w,
		print ""

main()

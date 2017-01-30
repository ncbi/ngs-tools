#!/usr/bin/python
import sys
import math
import glob

def load_texts(filename):
	seqs = []
	f = open(filename)
	for line in f:
		line = line.rstrip()
		if line:
			seqs.append(line)
		
	return seqs

class Kmer:
	def __init__(self, text, weight):
		self.text = text
		self.weight = weight

MISSING_CHAR_KMER = Kmer(' ', 20)

class Compression:
	def __init__(self, left_part, kmer):
		self.left_part = left_part
		self.kmer = kmer
		self._weight = self.calc_weight()

	def calc_weight(self):
		w = self.kmer.weight
		if self.left_part != None:
			w += self.left_part.weight()

		return w

	def weight(self):
		return self._weight

	def kmers(self):
		left = []
		left_part = self.left_part
		
		while left_part != None:
#			left = self.left_part.kmers()
#			left = [left_part.kmer] + left
			left.append(left_part.kmer)
			left_part = left_part.left_part

		left.reverse()

		left.append(self.kmer)
		return left # + [self.kmer]

class Kmers:
	def __init__(self):
		self.kmers = dict()
 
	def get_kmer(self, text):
		dic = self.kmers.get(len(text), None)
		if dic == None:
			return None

		kmer = dic.get(text, None)
		if len(text) == 1 and kmer == None:
			return MISSING_CHAR_KMER

		return kmer

	def add(self, text, weight_bits):
		if not len(text) in self.kmers:
			self.kmers[len(text)] = dict()

		self.kmers[len(text)][text] = Kmer(text, weight_bits)
#		print text, weight_bits

	def max_kmer_len(self):
		return max([length for length in self.kmers])

def load_kmers(file_mask):
	kmers = Kmers()
	for filename in glob.glob(file_mask):
		load_kmers_from_file(kmers, filename)

	return kmers

def load_kmers_from_file(kmers, filename):
	lines = []
	
	f = open(filename)
	for line in f:
		line = line.split('|')
		lines.append((line[0], int(line[1])))

todo: revert - calc sum weight for all lines, not for file
	sum_weight = sum([count for (kmer, count) in lines]) # todo: fix bug. why we calculate total sum of all kmers ????

	for kmer, count in lines:
#		print '|' + kmer + '|'
		kmers.add(kmer, bits_for(count, sum_weight))

def combine(compression, kmer):
	if kmer == None:
		return None
	return Compression(compression, kmer)

def get_kmer(kmers, text, kmer_len, _from):
	s = text[_from : _from + kmer_len]
	res = kmers.get_kmer(s)

#	print "getting kmer|" + s + '|',
#	if res != None:
#		print "ok"
#	else:
#		print "None"

	return res

LONG_KMER_LEN = 6

def print_compression(compressed):
	weight = 0
	kmers = compressed.kmers()
	for kmer in kmers:
		print '|' + kmer.text + '(' + str(round(kmer.weight, 2)) + ')' + '|',
#		print '|' + kmer.text + '|',
		weight += kmer.weight

	print " = ", weight

#	text_len = 0
#	long_mers_text_len = 0
#	for kmer in kmers:
#		text_len += len(kmer.text)
#		if len(kmer.text) >= LONG_KMER_LEN:
#			print kmer.text,
#			long_mers_text_len += len(kmer.text)

#	s = " = " + str(long_mers_text_len) + " of " + str(text_len) + " = " + str( round(100.0*long_mers_text_len/text_len, 2) ) + "%"
#	print s
#	sys.stderr.write(s + '\n')

def compress(text, kmers):
	compression = dict()

#	print "max kmer len is", kmers.max_kmer_len()
	for to_inclusive in range(len(text)):
		best_compression = None
#		if to_inclusive % 1024 == 0:
#			sys.stderr.write('.')
#		print "----------------------", to_inclusive
			
		for kmer_len in range(1, kmers.max_kmer_len() + 1):
			left_part_to_inclusive = to_inclusive - kmer_len
#			print "left part to", left_part_to_inclusive
			left_part = None
			if left_part_to_inclusive < -1:
				break

			if left_part_to_inclusive >= 0:
				left_part = compression[left_part_to_inclusive]

			compressed = combine(left_part, get_kmer(kmers, text, kmer_len, left_part_to_inclusive + 1))
#			if compressed != None:
#				print_compression(compressed)
			if compressed != None and (best_compression == None or compressed.weight() < best_compression.weight()):
				best_compression = compressed

		if best_compression == None:
			print "bad letter:", to_inclusive, text[to_inclusive]
			raise Exception("best_compression == None")

		compression[to_inclusive] = best_compression
#		print "best compression for", to_inclusive, "set"
#		print_compression(best_compression)

	return compression[len(text) - 1]

def calc_kmers_freq(real_kmers):
	freq = dict()
	for kmer in real_kmers:
		freq[kmer.text] = freq.get(kmer.text, 0) + 1

	return freq

def dump_frequency(real_kmers, loaded_kmers, min_coverage):
	kmers_freq = calc_kmers_freq(real_kmers)
	total_words = 0
	for length in loaded_kmers.kmers:
#		print length
		filename = "stat" + str(length) + ".txt"
		f = open(filename, 'w')
		f.flush()

		min_cov = min_coverage
		if length <= 2:
			min_cov = 0

		words_count = dump_frequency_of_len(f, kmers_freq, loaded_kmers.kmers[length], min_cov)
		if length >= LONG_KMER_LEN:
			total_words += words_count

	sys.stderr.write("total long words: " + str(total_words) + '\n')

def bits_for(count, sum_weight):
	MANY_BITS = 20
	if count == 0:
		return MANY_BITS
	part = 1.0*count/sum_weight
	return math.log(1.0/part, 2)
#	print kmer + '|' + str(count) + '|' + str( round(bits, 4) )


def dump_frequency_of_len(f, kmers_freq, kmers_of_len, min_coverage):
#	sum_freq = 0
#	for kmer_text in kmers_of_len:
#		sum_freq += kmers_freq.get(kmer_text, 0)

#	f.write(str(sum_freq) + '\n')

	lines = []
#	sys.stderr.write("dump next freq of\n")

	for kmer_text in kmers_of_len:
		times = kmers_freq.get(kmer_text, 0)
		if times >= min_coverage:
			lines.append((kmer_text + '|' + str(times), times))
#			lines.append((kmer_text + '|' + str(times) + '|' + str(round(bits_for(times, sum_freq), 4)), times))

	lines.sort(key = lambda( text, count) : -count)

	for text, count in lines:
		f.write(text + '\n')

	return len(lines)
	
def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 4:
		print "need <text file> <kmers path mask> <min dump coverage>"
		return

	texts = load_texts(sys.argv[1])
	kmers = load_kmers(sys.argv[2])

	results = []
	for text in texts:
		compressed = compress(text, kmers)
		print_compression(compressed)
		results.append(compressed)
		sys.stderr.write('.')
		
#	sys.stderr.write("dump freq\n")
#	dump_frequency(merge_lists([compressed.kmers() for compressed in results], kmers, int(sys.argv[3]))
	
main()

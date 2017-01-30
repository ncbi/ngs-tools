#!/opt/python-2.7/bin/python
#!/usr/bin/python
import os
import sys
import perf
import string
import shell

KMER_LEN = 32
KMER_CONTAMINATION_LEN = 16
MIN_COVERAGE = 2
CONTAMINATION_COVERAGE_MULT = 4
MIN_CONTIG_LEN = 200
MAX_CONTIGS = 1000
ABSOLUTE_MIN_CONTIG_LEN = 200
CONTAMINATION_COVERAGE_MEDIAN_FILTER_RADIUS = 5
#MIN_CONSENSUS_PERCENT = 90

def get_all_variants(seq):
	seq_compl = seq.translate(string.maketrans("ACTG", "TGAC"))
	return [seq, seq[::-1], seq_compl, seq_compl[::-1]]

def remove_all_variants(kmers, x):
	x = get_all_variants(x)
	for t in x:
		if t in kmers:
			del kmers[t]

def add_all_variants(kmers, x):
# todo: register only ACTG
#	all_variants = get_all_variants(x)
#	for x in all_variants:
	kmers[x] = kmers.get(x, 0) + 1
#		kmers[x[0]] = kmers.get(x[0], 0) + 1
	return 1

def update_dict_with_read(kmers, line, kmer_len):
	count = 0
	for i in range(len(line) - kmer_len + 1):
		kmer = line[i:i + kmer_len]
#		if len(kmer) != kmer_len:
#			raise Exception("len(kmer) != kmer_len")
		count += add_all_variants(kmers, kmer)

	return count

#@perf.timed		
def build_kmer_dict(filename, kmer_len, kmer2_len):
	kmers = dict()
	kmers2 = dict()
	f = open(filename)
	line_count = 0
	total_weight = 0
	for line in f:
		if line and line[0] != '>':
			line = line[:-1] # w/o \n at the end
			update_dict_with_read(kmers, line, kmer_len)
			total_weight += update_dict_with_read(kmers2, line, kmer2_len)
			line_count += 1
#			if line_count % 1000000 == 0:
#				print line_count
#				sys.stdout.flush()

	return kmers, kmers2, total_weight

def any_of(kmers):
	for x in kmers:
		return x

	return None

#@perf.timed
#def get_best_kmer(kmers):
#	if not kmers:
#		return None, None

	# todo: optimize
#	x = any_of(kmers)
#	return x, kmers[x]
			
#	best_score = 0
#	best_kmer = None
#	for x in kmers:
#		if kmers[x] > best_score:
#			best_score = kmers[x]
#			best_kmer = x

#	return best_kmer, best_score
		
#@perf.timed
def get_next_contig(kmers, kmer_len, seq):
#	seq, seq_cov = get_best_kmer(kmers)
#	if not seq or seq_cov < MIN_COVERAGE:
#		return None, None
	seq_cov = kmers[seq]

	remove_all_variants(kmers, seq)
	cov = [(seq_cov, seq_cov)] * len(seq)

	was_reversed = False

	while True:
		next_letter, letter_cov, letter_cov_sum = choose_next_letter(seq, kmers, kmer_len)
		if not next_letter: # or unexpected_coverage_behaviour(last some + letter_cov ...): todo: implement
			if not was_reversed:
				seq = seq[::-1]
				cov = cov[::-1]
				was_reversed = True
			else:
#				if was_reversed:
				seq = seq[::-1]
				cov = cov[::-1]
				return seq, cov
		else:
			seq += next_letter
			cov.append((letter_cov, letter_cov_sum))

def get_kmer_candidates(seq):
	return [seq + 'A', seq + 'C', seq + 'T', seq + 'G']

def choose_best_kmer_candidate(candidates):
	candidates.sort(key = lambda (seq, cov): - cov)
	sum_coverage = sum([cov for (seq, cov) in candidates])
	top_2 = candidates[0][1] + candidates[1][1]
	if sum_coverage >= MIN_COVERAGE: # and expected_coverage_change(prev_coverage, sum_coverage): # has_consensus(top_2, sum_scores) and :
		return candidates[0]
	return None, None

def coverage_of(kmers, seq):
	seqs = get_all_variants(seq)
	return sum([kmers.get(s, 0) for s in seqs])

#def get_seq_snps(seq, _from, _to):
#	seqs = []
#	letters = ['A', 'C', 'T', 'G']
#	for pos in range(_from, _to):
#		for l in letters:
#			if l != seq[pos]:
#				seqs.append(seq[:pos] + l + seq[pos + 1:])
#
#	return seqs

def coverage_of_snp(kmers, seq):
	cov = coverage_of(kmers, seq)
	direct_snps = get_seq_snps(seq, 0, len(seq) - 1)
	for snp_seq in direct_snps:
		cov = max(cov, coverage_of(kmers, snp_seq)) #kmers.get(snp_seq, 0)

	return cov

#	back_snps = get_seq_snps(seq[::-1], 0, len(seq) - 1)
#	for snp_seq in back_snps:
#		cov = max(cov, kmers.get(snp_seq, 0)
	

def choose_next_letter(seq, kmers, kmer_len):
	seq = seq[ - (kmer_len - 1) : ]
	candidates = get_kmer_candidates(seq)
	seq_cov = [(seq, coverage_of(kmers, seq)) for seq in candidates]

	c, letter_cov = choose_best_kmer_candidate(seq_cov)
	if not c:
		return None, None, None

	for x in candidates:
		remove_all_variants(kmers, x)
	
	return c[-1], letter_cov, sum([cov for (seq, cov) in seq_cov]) # todo: remove sum

#@perf.timed
#def clear_kmers(kmers, min_coverage):
#	to_remove = [] # todo: can we do it in place?
#	for x in kmers:
#		if kmers[x] < min_coverage:
#			to_remove.append(x)

#	for x in to_remove:
#		del kmers[x]

#def get_cov_delta(cov):
#	cov_delta = []
#	for i in range(1, len(cov)):
#		prev = cov[i-1]
#		current = cov[i]
#		cov_delta.append(100*current/prev)
#
#	return cov_delta

def choose_begin(begins, kmers, begin_i):
	if not kmers:
		return None

	begin_i += 1
	while begin_i < len(begins):
		seq, cov = begins[begin_i]
		if seq in kmers:
			return begin_i

		begin_i += 1

	return None

def calculate_coverage(contig, kmers, kmer_len):
	cov = []
	for i in range(len(contig) - kmer_len + 1):
		s = contig[i: i + kmer_len]
		seqs = get_all_variants(s)
		cov.append(sum([kmers.get(seq, 0) for seq in seqs]))

	return cov

def median_of(_v):
	v = list(_v)
	v.sort()
	return v[len(v)/2]

def split_sequence(seq, contigs):
	begins = [pos for c, pos in contigs] + [len(seq)]
	ends = [0] + [pos + len(c) for c, pos in contigs]
	result = [c for c, pos in contigs]

	for i in range(len(contigs) + 1):
		beg = ends[i]
		end = begins[i]

		if beg > end:
#			print seq
#			print contigs
			raise Exception("split_sequence beg > end")

		if beg == end:
			continue

		result.append(seq[beg:end])

	return result

class ContigCleaner:
	def __init__(self, contig, cont_coverage, kmer_len):
		self.contig = contig
		self.cont_coverage = cont_coverage
		self.avg_cov = median_of(cont_coverage)
		self.kmer_len = kmer_len
#		print "median is", self.avg_cov
		self.result = []
		self.contaminations = []

	def clean(self):
		self.good = []
		self.bad = []

		for i in range(len(self.contig)):
			if self.cont_cov(i) < self.avg_cov * CONTAMINATION_COVERAGE_MULT: # todo: smarter clustering - we are missing some repeats now
				self.belong_to_good(i)
			else:
				self.belong_to_bad(i)

		self.flush_good(self.good)
		self.flush_bad(self.bad)
		return split_sequence(self.contig, self.result)

	def belong_to_good(self, i):
		self.bad = self.flush_bad(self.bad)
		if not self.good:
			self.good_start_pos = i
		self.good.append(self.contig[i])
		
	def belong_to_bad(self, i):
		self.good = self.flush_good(self.good)
		self.bad.append(self.contig[i])

	def cont_cov(self, i):
		if i < len(self.cont_coverage):
			return self.cont_coverage[i]

		return self.cont_coverage[-1]

	def flush_bad(self, seq):
		if seq:
			self.contaminations.append(seq)
		return []

	def flush_good(self, seq):
		if seq:
			seq = ''.join(seq)
			if self.contaminations:
				seq = seq[self.kmer_len:]

			self.result.append((seq, self.good_start_pos))
				
		return []

def median_filter(v, radius):
	result = []
	for pos in range(radius/2):
		result.append(v[pos])

	for pos in range(len(v) - radius):
		rs = []
		for dx in range(radius):
			rs.append(v[pos + dx])
		x = median_of(rs)
		
		result.append(x)

	while len(result) < len(v):
		result.append(v[len(result)])

	return result

def build_contigs(fasta_filename):
	kmers, contaminated_kmers, total_weight = build_kmer_dict(fasta_filename, KMER_LEN, KMER_CONTAMINATION_LEN)
#	print "total_weight", total_weight
#	cov = [kmers[x] for x in kmers]
#	cov.sort()

#	sys.stdout.flush()
	begins = [(seq, kmers[seq])	for seq in kmers if kmers[seq] >= MIN_COVERAGE]
	begins.sort(key = lambda (seq, cov): - cov)
	begin_i = -1
#	print "size of begins", len(begins)
#	print begins[0][0], begins[0][1]

	result = []

	while True:
		begin_i = choose_begin(begins, kmers, begin_i)
		if begin_i == None:
			break

		contig, cov = get_next_contig(kmers, KMER_LEN, begins[begin_i][0]) # todo: coverage don't need I think
		if not contig:
			break

		if len(contig) < MIN_CONTIG_LEN:
			continue

		cleaner = ContigCleaner(contig, median_filter(calculate_coverage(contig, contaminated_kmers, KMER_CONTAMINATION_LEN), CONTAMINATION_COVERAGE_MEDIAN_FILTER_RADIUS), KMER_CONTAMINATION_LEN)
# todo: remove next 2 lines
		cov = calculate_coverage(contig, contaminated_kmers, KMER_CONTAMINATION_LEN)
		cov = median_filter(cov, CONTAMINATION_COVERAGE_MEDIAN_FILTER_RADIUS)

		contigs = cleaner.clean() 
		for c in contigs:
			if len(c) < ABSOLUTE_MIN_CONTIG_LEN:
				continue
			cov = calculate_coverage(c, contaminated_kmers, KMER_CONTAMINATION_LEN)
			average_cov = sum(cov)/len(cov)
			data_percent = 100.0*sum(cov)/total_weight
			data_percent_str = '%.2f' % data_percent 
			desc = data_percent_str + "%_cov_" + str(average_cov) + "_len_" + str(len(c))
#			print c
			result.append((data_percent, c, desc))
#			if len(result) >= MAX_CONTIGS:
#				return result
		
#		sys.stdout.flush()

	return result


#@perf.timed
def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <fasta file>" # [-align]"
		return

	filename = sys.argv[1]
	delete_file = False
	if not '.' in filename:
		acc = sys.argv[1]
		shell.execute_local("fastq-dump " + acc + " --split-spot --skip-technical --fasta 0", False)
		filename = "./" + acc + ".fasta"
		delete_file = True

	result = build_contigs(filename)
	result.sort(key = lambda (data_percent, c, desc): - data_percent)
	result = result[:MAX_CONTIGS]
	index = 0	
	for data_percent, c, desc in result:
		print ">" + str(index) + '_' + desc
		print c
		index += 1

	if delete_file:
		os.remove(filename)

main()

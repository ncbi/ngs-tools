#!/opt/python-2.7/bin/python
import contig_builder
import perf

def equal(a, b):
	if a!=b:
		print str(a) + " != " + str(b)
	else:
		print "ok"

@perf.timed
def main():
#	kmers = contig_builder.build_kmer_dict("./spots/SRR1553418.fasta")
#	print len(kmers)
#	return
#	print contig_builder.split_sequence("GCTACGGACCTCCACCAGAGTTTCCTCTGGCTTCGCCCTGCCCAGGCATAGTTCCTGTCTCTTATACACATCTGACGCTGCCGACGAGCGATCTAGTGT", ['GCTACGGACCTCCACCAGAGTTTCCTCTGGCTTCGCCCTGCCCAGGCATAGTT', 'T'])
	print contig_builder.split_sequence("firstsecondthird", [("firstsecondthird", 0)])
	print contig_builder.split_sequence("firstsecondthird", [("first", 0)])
	print contig_builder.split_sequence("firstsecondthird", [("second", 5)])
	print contig_builder.split_sequence("firstsecondthirdfourthfifth", [("first", 0), ("third", 11), ("fourth", 16)])
	print contig_builder.split_sequence("firstsecondthirdfourthfifth", [("first", 0), ("third", 11), ("fifth", 22)])
	return


	print contig_builder.median_filter([1, 2, 3, -40, 5, 6, 7, 8, 9, 10, 40, 11, 244, 5, 2, -10, 4, 5], 5)
	return

	equal(contig_builder.get_all_variants("ACAG"), ["ACAG", "GACA", "TGTC", "CTGT"])
#	equal(contig_builder.get_all_variants("ACNAG"), [])

	kmers = dict()
	contig_builder.add_all_variants(kmers, "ACAG")
	equal(len(kmers), 1)
	contig_builder.add_all_variants(kmers, "GACA")
	equal(len(kmers), 2)
	contig_builder.remove_all_variants(kmers, "ACAG")
	equal(len(kmers), 0)
	return
	


	seq = "ACTGAA"
	print contig_builder.get_seq_snps(seq, 0, len(seq) - 1)
	return

	equal(contig_builder.choose_best_kmer_candidate([('a', 1), ('c', 2), ('t', 3), ('g', 0)]), ('t', 3))
	equal(contig_builder.choose_best_kmer_candidate([('a', 0), ('c', 0), ('t', 0), ('g', 0)]), (None, None))
	equal(contig_builder.choose_best_kmer_candidate([('a', 3), ('c', 0), ('t', 0), ('g', 0)]), (None, None))
	equal(contig_builder.choose_best_kmer_candidate([('a', 3), ('c', 1), ('t', 0), ('g', 0)]), ('a', 3))
	equal(contig_builder.choose_best_kmer_candidate([('a', 0), ('c', 2), ('t', 3), ('g', 0)]), ('t', 3))
	equal(contig_builder.choose_best_kmer_candidate([('a', 1), ('c', 2), ('t', 4), ('g', 3)]), (None, None))
	equal(contig_builder.choose_best_kmer_candidate([('a', 2), ('c', 20), ('t', 3), ('g', 0)]), ('c', 20))
#	print contig_builder.choose_best_kmer_candidate([('a', 0), ('c', 0), ('t', 0), ('g', 0)])
	return


	kmers = dict()
	contig_builder.update_dict_with_read(kmers, "CTATGGA", 4)
	contig_builder.update_dict_with_read(kmers, "TATGGAT", 4)
	print kmers
	print contig_builder.get_next_contig(kmers, 4)
	return

	print contig_builder.get_best_kmer(kmers)
	print contig_builder.choose_next_letter("GGTA", kmers, 4)
	return
	
	

	
	

main()
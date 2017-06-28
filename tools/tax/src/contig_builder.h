#ifndef CONTIG_BUILDER_H_INCLUDED
#define CONTIG_BUILDER_H_INCLUDED

#include <iomanip>
#include <numeric>
#include <sstream>
#include <math.h>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <chrono>
#include <thread>
#include "omp_adapter.h"

#include "hash.h"

struct ContigBuilder
{
	static const int MIN_COVERAGE = 2;
	static const int LAST_LETTERS_COUNT = 4;

	template <class KmerMap>
	static char choose_next_letter(KmerMap &kmers, typename KmerMap::hash_t *_hash, int min_coverage = MIN_COVERAGE)
	{
		auto hash = *_hash;
		const char LAST_LETTERS[LAST_LETTERS_COUNT] = {'A', 'C', 'T', 'G'};
		typename KmerMap::hash_t hashes[LAST_LETTERS_COUNT];
		unsigned int cov[LAST_LETTERS_COUNT]; 
	
        #pragma omp parallel for num_threads(LAST_LETTERS_COUNT)
		for (int i = 0; i < LAST_LETTERS_COUNT; i++)
		{
			hashes[i] = Hash<typename KmerMap::hash_t>::hash_next(LAST_LETTERS[i], hash, kmers.kmer_len);
			cov[i] = kmers.coverage_of(hashes[i]); 
		}

		int best_letter_index = 0;
		unsigned int best_letter_cov = cov[best_letter_index];
		unsigned int sum_cov = best_letter_cov;

		for (int i = 1; i < LAST_LETTERS_COUNT; i++)
		{
			unsigned int current_cov = cov[i];
			sum_cov += current_cov;
			if (current_cov > best_letter_cov)
			{
				best_letter_cov = current_cov;
				best_letter_index = i;
			}
		}

		if (best_letter_cov < min_coverage)
			return 0;

		*_hash = hashes[best_letter_index];
		return LAST_LETTERS[best_letter_index];
	}

	template <class KmerMap>
	static std::string get_next_contig(KmerMap &kmers, typename KmerMap::hash_t start_from, int min_coverage = MIN_COVERAGE)
	{
		std::string seq = Hash<typename KmerMap::hash_t>::str_from_hash(start_from, kmers.kmer_len);
		auto hash = start_from;
		kmers.remove(hash);

		bool was_reversed = false;

		while (true)
		{
			char next_letter = choose_next_letter<KmerMap>(kmers, &hash, min_coverage);
			if (!next_letter)
			{
				if (!was_reversed)
				{
                    seq_transform_actg::to_rev_complement(seq);
					hash = seq_transform<typename KmerMap::hash_t>::to_rev_complement(start_from, kmers.kmer_len);
					was_reversed = true;
				}
				else
				{
                    seq_transform_actg::to_rev_complement(seq);
					return seq;
				}
			}
			else
			{
				kmers.remove(hash);
				seq += next_letter;
			}
		}
	}

};

#endif

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
#include <omp.h>

#include "hash.h"
//#include "contig_cleaner.h"
//#include "kmer.h"
//#include "coverage.h"
//#include "contig.h"
//#include "run_reader.h"

struct ContigBuilder
{
	static bool almost_equal(const std::string &best, int best_from, const std::string &second, int second_from, int size, int max_errors = 3)
	{
		for (int i=0; i<size; i++)
			if (best[i + best_from] != second[i + second_from])
			{
				max_errors--;
				if (max_errors < 0)
					return false;
			}

		return true;
	}


	static bool small_difference_at_beg(const std::string &best, const std::string &second)
	{
		return 
			almost_equal(best, 0, second, 0, best.size()/2) ||
			almost_equal(best, 1, second, 0, best.size()/2) ||
			almost_equal(best, 0, second, 1, best.size()/2) ||
			almost_equal(best, 2, second, 0, best.size()/2) ||
			almost_equal(best, 0, second, 2, best.size()/2);
	}

	static bool small_difference_at_end(const std::string &best, const std::string &second)
	{
		return 
			almost_equal(best, best.size()/2, second, second.size()/2, best.size()/2) ||
			almost_equal(best, best.size()/2 - 1, second, second.size()/2, best.size()/2) ||
			almost_equal(best, best.size()/2, second, second.size()/2 - 1, best.size()/2) ||
			almost_equal(best, best.size()/2 - 2, second, second.size()/2, best.size()/2) ||
			almost_equal(best, best.size()/2, second, second.size()/2 - 2, best.size()/2);
	}

	static const int MIN_COVERAGE = 2;

	template <class KmerMap>
	static bool cov_consensus(KmerMap &kmers, int best_letter_index, unsigned int sum_cov, unsigned int *cov, typename KmerMap::hash_t *hashes, int min_coverage)
	{
		int best_cov = cov[best_letter_index];
		if (sum_cov - best_cov < min_coverage)
			return true;

		int second_cov = 0;
		int second_cov_index = -1;
		for (int i=0; i<LAST_LETTERS_COUNT; i++)
		{
			unsigned int current_cov = cov[i];
			if (i != best_letter_index && current_cov > second_cov)
			{
				second_cov = current_cov;
				second_cov_index = i;
			}
		}

		if (second_cov < min_coverage)
			return true;

		int test_len = kmers.kmer_len; // todo: tune

		typename KmerMap::hash_t second_hash = hashes[second_cov_index];
		for (int i=0; i<test_len; i++)
			if (!choose_next_letter(kmers, &second_hash, std::max(second_cov/2, min_coverage), false))
				return true;

		typename KmerMap::hash_t best_hash = hashes[best_letter_index];
		for (int i=0; i<test_len; i++)
			if (!choose_next_letter(kmers, &best_hash, std::max(best_cov/2, min_coverage), false))
				return true;

		if (second_hash == best_hash) // does not matter which letter to choose
			return true;

		// todo: optimize
		std::string best = Hash<typename KmerMap::hash_t>::str_from_hash(best_hash, kmers.kmer_len);
		std::string second = Hash<typename KmerMap::hash_t>::str_from_hash(second_hash, kmers.kmer_len);

//		if (small_difference_at_end(best, second))
		if (small_difference_at_beg(best, second))
			return true;

		//std::cout << "cov: " << best_cov << "/" << sum_cov << std::endl;
		//std::cout << "best   " << best << std::endl;
		//std::cout << "second " << second << std::endl;
		return false;
	}

	static const int LAST_LETTERS_COUNT = 4;
	template <class KmerMap>
	static char choose_next_letter(KmerMap &kmers, typename KmerMap::hash_t *_hash, int min_coverage = MIN_COVERAGE, bool check_consensus = true)
	{
		auto hash = *_hash;
		const char LAST_LETTERS[LAST_LETTERS_COUNT] = {'A', 'C', 'T', 'G'};
		typename KmerMap::hash_t hashes[LAST_LETTERS_COUNT];
		unsigned int cov[LAST_LETTERS_COUNT]; 
	
        #pragma omp parallel for num_threads(LAST_LETTERS_COUNT)
		for (int i=0; i<LAST_LETTERS_COUNT; i++)
		{
			hashes[i] = Hash<typename KmerMap::hash_t>::hash_next(LAST_LETTERS[i], hash, kmers.kmer_len);
			cov[i] = kmers.coverage_of(hashes[i]); // todo: insert barrier next ?
//			remove_all_variants(kmers, hashes[i]); //seq + LAST_LETTERS[i]);
		}

		int best_letter_index = 0;
		unsigned int best_letter_cov = cov[best_letter_index];
		unsigned int sum_cov = best_letter_cov;

		for (int i=1; i<LAST_LETTERS_COUNT; i++)
		{
			unsigned int current_cov = cov[i];
			sum_cov += current_cov;
			if (current_cov > best_letter_cov)
			{
				best_letter_cov = current_cov;
				best_letter_index = i;
			}
		}

//		std::cout << "sum cov: " << sum_cov << " of min cov " << min_coverage << std::endl;
//		std::cout << "cov: " << best_letter_cov << "/" << sum_cov << std::endl;
		if (best_letter_cov < min_coverage || (check_consensus && !cov_consensus(kmers, best_letter_index, sum_cov, cov, hashes, min_coverage)))
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
			char next_letter = choose_next_letter<KmerMap>(kmers, &hash, min_coverage, false);
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

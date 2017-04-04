#ifndef KMER_MAP_H_INCLUDED
#define KMER_MAP_H_INCLUDED

#include <unordered_map>
#include <map>
#include <mutex>
#include <iostream>
#include "seq_transform.h"
#include "hash.h"
#include <omp.h>

template <class _hash_t, int _kmer_len, int _count_buckets>
struct KmerMap 
{
	typedef _hash_t hash_t;
	static const int COUNT_BUCKETS = _count_buckets;
	static const int kmer_len = _kmer_len;

	std::mutex kmer_mutex;

	struct Count
	{
		unsigned int reverse : 1;
		unsigned int complement : 1;
		unsigned int deleted : 1;
		unsigned int count : 29;
		static const unsigned int MAX_COUNT = (1 << 29) - 1;

		Count(int count = 0) : count(count), reverse(0), complement(0), deleted(0){} // todo: check performance
	};

	typedef std::unordered_map<hash_t, Count> CountMap;
//	typedef std::map<hash_t, Count> CountMap;
	std::vector< CountMap > count;
	std::vector< std::mutex > bucket_mutex;
	std::vector<long long unsigned int> bucket_weight;
	std::vector<size_t> bucket_frequent;

	KmerMap() : count(COUNT_BUCKETS), bucket_mutex(COUNT_BUCKETS), bucket_weight(COUNT_BUCKETS), bucket_frequent(COUNT_BUCKETS)
	{
		static_assert(sizeof(Count) == sizeof(int), "sizeof(Count) == sizeof(int)");
	}

	void add(hash_t hash)
	{
		bool complement = false, reverse = false;
		hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len, &complement, &reverse);

		auto bucket = get_count_bucket(hash);
		bucket_mutex[bucket].lock(); // todo: try atomic?
		auto &c = count[bucket][hash];

		if (c.count == 0)
		{
			c.complement = complement;
			c.reverse = reverse;
		}

		if (c.count == 1)
			bucket_frequent[bucket]++;

		if (c.count < Count::MAX_COUNT)
			c.count++;

		bucket_weight[bucket]++;
		bucket_mutex[bucket].unlock();
	}

	void reserve(size_t size, float max_load_factor = 4.0f) // todo: tune max load factor
	{
//		const float MAX_LOAD_FACTOR = 4.0f;
//		std::cerr << "reserve " << size/count.size() << " for each" << std::endl;
		for (auto &c : count)
		{
			c.max_load_factor(max_load_factor);
			c.reserve(size/count.size());
		}
	}

	static unsigned int get_count_bucket(hash_t hash)
	{
		auto bucket = (hash >> 2) % COUNT_BUCKETS; // different last letter sequence should have the same bucket
		return bucket;
	}

	unsigned int get(hash_t hash) const
	{
		auto c = get_full(hash);
		return c.deleted ? 0 : c.count;
	}

	Count get_full(hash_t hash) const
	{
		auto &bucket = count[get_count_bucket(hash)];
		auto it = bucket.find(hash);
		if (it == bucket.end()) // || it->second.deleted)
			return Count();

		return it->second; //.count;
	}

	void remove(hash_t hash)
	{
		bool complement = false, reverse = false;
		hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len, &complement, &reverse);

		auto &bucket = count[get_count_bucket(hash)];
		auto it = bucket.find(hash);
		if (it == bucket.end() || it->second.deleted)
			return;

		it->second.deleted = 1;
	}

	void restore(hash_t hash)
	{
		bool complement = false, reverse = false;
		hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len, &complement, &reverse);

		auto &bucket = count[get_count_bucket(hash)];
		auto it = bucket.find(hash);
		if (it == bucket.end() || !it->second.deleted)
			return;

		it->second.deleted = 0;
	}

	unsigned int coverage_of(hash_t hash) const
	{
		return get(seq_transform<hash_t>::min_hash_variant(hash, kmer_len));
	}

	unsigned int coverage_of_no_deleted_check(hash_t hash) const
	{
		return get_full(seq_transform<hash_t>::min_hash_variant(hash, kmer_len)).count;
	}

	void get_original_compl_rev(hash_t hash, bool *orig_complement, bool *orig_reverse) const
	{
		bool complement = false, reverse = false;
		hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len, &complement, &reverse);
		auto c = get_full(hash);
//		if (!c.count)
	//		throw std::runtime_error("originally_reversed test for non-existing");

		*orig_complement = complement != c.complement; 
		*orig_reverse = reverse != c.reverse; 
	}

	bool originally_complement(hash_t hash) const
	{
		bool complement = false, reverse = false;
		get_original_compl_rev(hash, &complement, &reverse);
		return complement;
	}

	bool originally_reverse(hash_t hash) const
	{
		bool complement = false, reverse = false;
		get_original_compl_rev(hash, &complement, &reverse);
		return reverse;
	}

	long long unsigned int total_weight() const
	{
		long long unsigned int sum = 0;
		for (auto &w : bucket_weight)
			sum += w;

		return sum;
	}

	size_t size() const
	{
		size_t sum = 0;
		for (auto &bucket : count)
			sum += bucket.size();
		
		return sum;
	}
			
	void optimize(int min_count = 2) // todo: multitheading
	{
		const float MAX_OPTIMIZED_LOAD_FACTOR = 4.0f;
		#pragma omp parallel for
		for (int bucket_i = 0; bucket_i < COUNT_BUCKETS; bucket_i++)
		{
			CountMap new_bucket;
			new_bucket.max_load_factor(MAX_OPTIMIZED_LOAD_FACTOR);
			new_bucket.reserve(bucket_frequent[bucket_i]);
//			std::cerr << "reserving " << bucket_frequent[bucket_i] << std::endl;
			auto &bucket = count[bucket_i];
			bucket_weight[bucket_i] = 0;

			for (auto c = bucket.begin(); c != bucket.end(); c++)
				if (c->second.count >= min_count)
				{
					new_bucket.insert(*c);
					bucket_weight[bucket_i] += c->second.count;
				}

//			std::cerr << "result size " << new_bucket.size() << std::endl;
			bucket.swap(new_bucket);
		}
	}

	template <class Lambda>
	void for_every_kmer_do(Lambda &&lambda) const // todo: decide what to do with deleted
	{
		for (auto &bucket : count)
		for (auto &map_element: bucket)
			if (!map_element.second.deleted)
				lambda(map_element.first, map_element.second.count);
	}

};

typedef KmerMap<__uint128_t, 64, 1024> KmerMap64;
//typedef KmerMap<__uint64_t, 32, 1024> KmerMap32;
typedef KmerMap<unsigned int, 16, 1024> KmerMap16;

#endif

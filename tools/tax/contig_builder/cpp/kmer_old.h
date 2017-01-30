#ifndef KMER_H_INCLUDED
#define KMER_H_INCLUDED

#include <string>
#include <map>
#include <unordered_map>
#include <list>
#include <algorithm>
#include <iostream>
#include <mutex>
#include <thread>
#include <chrono>
#include "seq_transform.h"
#include "hash.h"

struct Kmers;
static long long int load_reads_thread(Kmers &kmers, int from, int step);

static const int LOADING_THREADS = 8; // todo: move inside class - only for perf bug
struct Kmers 
{
//	static const int COUNT_BUCKETS_BITS = 9;
	static const int COUNT_BUCKETS = 1 << 9; //1 << COUNT_BUCKETS_BITS;

	std::mutex kmer_mutex;

	int kmer_len;
	std::vector< std::unordered_map<hash_t, int> > count;
	std::vector< std::mutex > count_mutex;
	long long int total_weight;

	Kmers(int kmer_len) : 
		kmer_len(kmer_len), total_weight(0), count(COUNT_BUCKETS), count_mutex(COUNT_BUCKETS)
	{
	}

	void load(const std::vector<std::string> &reads)
	{
		auto before = std::chrono::high_resolution_clock::now();

		{
			size_t expected_kmers_count = get_expected_kmers_count(reads);
			for (auto &c : count)
				c.reserve(expected_kmers_count/COUNT_BUCKETS);
			std::cerr << "expected kmers: " << expected_kmers_count << std::endl;
		}

		std::list<std::thread> threads;
		for (int i=0; i<LOADING_THREADS; i++)
			threads.push_back(std::thread(&Kmers::load_reads_thread, this, i, LOADING_THREADS, reads));

		for (auto &t : threads)
			t.join();

		print_count_stat();
		std::cerr << "building kmers dict time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - before ).count() << std::endl;
	}

private:

	size_t get_expected_kmers_count(const std::vector<std::string> &reads)
	{
		size_t kmers_count = 0;
		size_t max_kmers = size_t(1) << 32;//Hash::hash_values(kmer_len); 

		for (auto &r : reads)
		{
			int kmers = r.length() - kmer_len + 1;
			kmers_count += kmers;
			if (kmers_count >= max_kmers)
				return max_kmers;
		}

		return kmers_count;
	}


	void load_reads_thread(int from, int step, const std::vector<std::string> &reads)
	{
		long long int total_weight = 0; // todo: don't need it at all - just count of add() calls
		for (unsigned int i = from; i<reads.size(); i+=step)
			total_weight += update_with_read(reads[i]);

		kmer_mutex.lock();
		this->total_weight += total_weight;
		kmer_mutex.unlock();
	}

public:

	int update_with_read(const std::string &s)
	{
		if (s.length() < kmer_len)
			return 0;

		size_t hash = Hash::hash_of(&s[0], kmer_len);
		int weight = 0;
		for (int i=0; i < s.size() - kmer_len + 1; i++, hash = Hash::hash_next(&s[i], hash, kmer_len)) // todo: <= and remove +1
		{
			add(hash);
			weight++; // actually don't need it in this case
		}

		return weight;
	}

	void add(hash_t hash)
	{
		auto bucket = get_count_bucket(hash);
		count_mutex[bucket].lock();
		count[bucket][hash]++;
		count_mutex[bucket].unlock();
	}

	unsigned int get_count_bucket(hash_t hash)
	{
//		auto shift_bits = Hash::hash_bits(kmer_len) - COUNT_BUCKETS_BITS;
		auto bucket = (hash >> 2) % COUNT_BUCKETS; // different last letter sequence should have the same bucket
		return bucket;
	}

	unsigned int get(const std::string &seq)
	{
		return get(&seq[0]);
	}

	unsigned int get(const char *s)
	{
		size_t hash = Hash::hash_of(&s[0], kmer_len);
		return get(hash);
	}

	unsigned int get(hash_t hash)
	{
		auto &bucket = count[get_count_bucket(hash)];
		auto it = bucket.find(hash);
		if (it == bucket.end() || erased(it->second))
			return 0;

		return it->second;
	}

	unsigned int get_erased(hash_t hash)
	{
		auto &bucket = count[get_count_bucket(hash)];
		auto it = bucket.find(hash);
		if (it == bucket.end() || !erased(it->second))
			return 0;

		return - it->second;
	}


	void erase(const std::string &seq)
	{
		return erase(&seq[0]);
	}

	void erase(const char *s)
	{
		erase(Hash::hash_of(s, kmer_len));
	}

	void erase(hash_t hash)
	{
		auto &bucket = count[get_count_bucket(hash)];
		auto it = bucket.find(hash);
		if (it == bucket.end())
			return;

		if (!erased(it->second))
			it->second = -it->second; // hack to erase and keep it at the same moment

//		count[get_count_bucket(hash)][hash] = 0; // todo: optimize
	}

	bool erased(int count) const
	{
		return count <= 0;
	}

	unsigned int coverage_of(const std::string &seq) // todo: optimize
	{
		return coverage_of(&seq[0]);
	}

	unsigned int coverage_of(const char *s)
	{
		return coverage_of(Hash::hash_of(s, kmer_len)); 
	}

	unsigned int coverage_of(hash_t hash) // todo: optimize
	{
//		return coverage_of(Hash::str_from_hash(hash, kmer_len));
		unsigned int sum = 0;
		seq_transform::for_every_hash_variant_do(hash, kmer_len, [&sum, this](hash_t h){ sum += get(h); }); //todo: make hash variants too
		return sum;
	}

	unsigned int erased_coverage_of(const char *s)
	{
		return erased_coverage_of(Hash::hash_of(s, kmer_len)); 
	}

	unsigned int erased_coverage_of(hash_t hash) // todo: optimize
	{
		unsigned int sum = 0;
		seq_transform::for_every_hash_variant_do(hash, kmer_len, [&sum, this](hash_t h){ sum += get_erased(h); }); //todo: make hash variants too
		return sum;
	}

	long long int current_total_weight() // todo: optimize
	{
		long long int weight = 0;
		for_every_kmer_do([&weight](hash_t hash, int count)
		{
			weight += count;
		});

		return weight;
	}

	void print_count_stat() // todo: optimize
	{
		int count1 = 0, count_many = 0;
		for_every_kmer_do([&count1, &count_many](hash_t hash, int count)
		{
			if (count == 1)
				count1++;
			else
				count_many++;
		});

		std::cerr << "count == 1 for " << count1 << std::endl;
		std::cerr << "count != 1 for " << count_many << std::endl;
	}


	template <class Lambda>
	void for_every_kmer_do(Lambda &&lambda)
	{
		for (auto &bucket : count)
		for (auto &map_element: bucket)
			if (!erased(map_element.second))
				lambda(map_element.first, map_element.second);
	}
};

static uint64_t hash_left(__uint128_t x)
{
	union 
	{
		__uint128_t x128;
		struct
		{
			uint64_t lo;
			uint64_t hi;
		} x64;
	} u;

	u.x128 = x;

//	std::cout << u.x64.hi << std::endl;
//	std::cout << u.x64.lo << std::endl;

	return u.x64.hi;
}

static uint32_t hash_left(uint64_t x)
{
	union 
	{
		uint64_t x64;
		struct
		{
			int lo;
			int hi;
		} x32;
	} u;

	u.x64 = x;

	return u.x32.hi;
}




#endif

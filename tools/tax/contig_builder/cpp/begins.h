#ifndef BEGINS_H_INCLUDED
#define BEGINS_H_INCLUDED

#include <vector>

template <class KmerMap>
struct Begins
{
	using hash_t = typename KmerMap::hash_t;
	struct Beg
	{
		hash_t hash;
		unsigned int count;
		Beg(hash_t hash, unsigned int count) : hash(hash), count(count){}
	};

	std::vector<Beg> begins;
	const KmerMap &kmers;
	typename std::vector<Beg>::const_iterator begin_it;

	Begins(const KmerMap &kmers, size_t min_coverage = 0) : kmers(kmers)
	{
		begins.reserve(kmers.size());
		kmers.for_every_kmer_do( [&] (hash_t hash, unsigned int count) // todo: can be multithreaded
			{
				if (count >= min_coverage)
					begins.push_back(Beg(hash, count));
			});

		std::sort(begins.begin(), begins.end(), [&kmers](const Beg & a, const Beg &b)
			{
//				auto a_count = kmers.get(a); // not coverage_of() - just get()
//				auto b_count = kmers.get(b);
				return b.count == a.count ? b.hash < a.hash : (b.count < a.count);
			});

		begin_it = begins.begin();
	};

	bool next(hash_t *hash_result)
	{
		while (begin_it != begins.end())
		{
			auto hash = begin_it->hash;
			if (kmers.get(hash) > 0) // not coverage - just get
			{
				*hash_result = hash;
				begin_it ++;
				return true;
			}

			begin_it ++;
		};

		return false;
	}

};


#endif

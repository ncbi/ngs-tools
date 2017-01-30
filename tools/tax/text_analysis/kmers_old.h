#ifndef KMERS_H_INCLUDED
#define KMERS_H_INCLUDED

#include <map>
#include <vector>
#include <string>
#include "stringn.h"

struct Kmer
{
	const char *seq;
	float bits;

	Kmer(const char *seq, float bits) : seq(seq), bits(bits){}

	bool empty() const 
	{
		return !seq && bits == 0; // *this != MISSING_KMER
	}
};

const static Kmer MISSING_KMER = Kmer(nullptr, 0);

// todo: make it non public
//bool operator <(const char *s1, const std::string &s2)
//{
//	return string_compare(s1, s2.c_str(), s2.size()) < 0;
//}


struct Kmers
{
//	struct Comparator
//	{
//		//bool operator()(char const *lhs, char const *rhs) const
//		//{
//		//	return std::strcmp(lhs, rhs) < 0;
//		//}
//		bool operator()(const std::string s1, const std::string &s2) const
//		{
////			return std::strcmp(lhs, rhs) < 0;
//			return s1 < s2;
//		}
//
//		bool operator ()(const char *s1, const std::string &s2)
//		{
//			return string_compare(s1, s2.c_str(), s2.size()) < 0;
//		}
//
//		bool operator ()(const std::string &s1, const char *s2)
//		{
//			return string_compare(s1.c_str(), s2, s1.size()) < 0;
//		}
//	};

	typedef std::map<std::string, float> KmersToBitsMap;

	std::vector<KmersToBitsMap> kmers;
	void add(const std::string &s, float bits)
	{
		if (s.length() > kmers.size())
			kmers.resize(s.length() + 1);

		auto &storage = kmers[s.length()];
		storage[s] = bits;
	}

	size_t max_kmer_len() const
	{
		return kmers.size() - 1;
	}

	Kmer get_kmer(const char *s, unsigned int kmer_len) const
	{
		if (kmer_len >= kmers.size())
			return MISSING_KMER;

		auto &storage = kmers[kmer_len];
		auto it = storage.find(std::string(s, s + kmer_len)); // todo: optimize
		if (it == storage.end())
			return MISSING_KMER;

		return Kmer(s, it->second);
	}
};

#endif
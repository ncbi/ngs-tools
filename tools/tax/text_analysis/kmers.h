#ifndef KMERS_H_INCLUDED
#define KMERS_H_INCLUDED

#include <map>
#include <vector>
#include <string>
#include "stringn.h"
#include "p_string.h"

struct Kmer
{
	const char *seq;
	float bits;
	int len;

	Kmer(const char *seq, float bits, int len) : seq(seq), bits(bits), len(len){}

	bool empty() const 
	{
		return !seq && !bits && !len; // *this != MISSING_KMER
	}
};

const static Kmer MISSING_KMER = Kmer(nullptr, 0, 0);
const static Kmer MISSING_ONE_MER = Kmer(nullptr, 20, 1);

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

	typedef std::map<p_string, float> KmersToBitsMap;

	std::vector<KmersToBitsMap> kmers;
	void add(const char *s, int len, float bits)
	{
//		std::cout << s << " " << len << " " << bits << std::endl;
		if (len >= kmers.size())
			kmers.resize(len + 1);

		auto &storage = kmers[len];
		storage[p_string(s, len)] = bits;
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
//		auto it = storage.find(std::string(s, s + kmer_len)); // todo: optimize
		auto it = storage.find(p_string(s, kmer_len)); // todo: optimize
		if (it == storage.end())
//			return kmer_len == 1 ? MISSING_ONE_MER : MISSING_KMER;
			return MISSING_KMER;

		return Kmer(it->first.s, it->second, kmer_len);
	}
};

#endif
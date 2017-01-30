#ifndef COVERAGE_H_INCLUDED
#define COVERAGE_H_INCLUDED

#include <string>
#include <vector>

template <class KmerMap>
struct Coverage : public std::vector<int>
{
	Coverage(const std::string &contig, const KmerMap &kmers) 
	{
		reserve(contig.length());
 		for (int i=0; i <= int(contig.length()) - kmers.kmer_len; i++)
		{
			auto hash = Hash<typename KmerMap::hash_t>::hash_of(&contig[i], kmers.kmer_len);
			push_back(kmers.coverage_of_no_deleted_check(hash));
		}
	}
};

//struct ErasedCoverage : public std::vector<int>
//{
//	ErasedCoverage(const std::string &contig, Kmers &kmers) 
//	{
//		reserve(contig.length());
// 		for (int i=0; i < int(contig.length()) - kmers.kmer_len + 1; i++)
//			push_back(kmers.erased_coverage_of(&contig[i]));
//	}
//};
//

#endif

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <iostream>
#include <chrono>
#include "omp_adapter.h"
//#include "hash_sorted_array.h"

#include "log.h"
#include "dbs.h"
#include "fasta.h"
#include "hash.h"
#include "config_build_centroids.h"
#include "seq_transform.h"
#include "file_list_loader.h"

typedef uint64_t hash_t;
typedef std::vector<hash_t> HashSortedArray;

bool has_hash(hash_t hash, HashSortedArray &hash_array)
{
    auto first = hash_array.begin();
    auto last = hash_array.end();
    first = std::lower_bound(first, last, hash);
    return !((first == last) || (hash < *first) );
}

struct KmerCount
{
	hash_t kmer = 0;
	int count = 0;
	
	KmerCount(hash_t kmer, int count) : kmer(kmer), count(count){}

	bool operator < (const KmerCount &x) const
	{
		return kmer < x.kmer;
	}
};

using namespace std;
using namespace std::chrono;

const string VERSION = "0.15";

struct Centroid
{
	vector<KmerCount> kmers;
	size_t seq_count = 0;

    Centroid(map<hash_t, int> &all_kmers, size_t seq_count) : seq_count(seq_count)
	{
		kmers.reserve(all_kmers.size());
		for (auto &x : all_kmers)
			kmers.push_back(KmerCount(x.first, x.second));
		std::sort(kmers.begin(), kmers.end());
	}
};

template <class CanUseKmer>
Centroid build_centroid(const string &filename, int kmer_len, CanUseKmer &&can_use_kmer)
{
    Fasta fasta(filename);
	map<hash_t, int> all_kmers;
	size_t seq_count = 0;

    string seq;
    while (fasta.get_next_sequence(seq))
    {
        string desc = fasta.sequence_description();

		seq_count++;
		set<hash_t> kmers;
        Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
        {
			auto kmer = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
			if (can_use_kmer(kmer))
				kmers.insert(kmer); // dont want to count duplicate kmers
            return true;
        });

		for (auto &kmer : kmers)
			all_kmers[kmer]++;
    }    

    return Centroid(all_kmers, seq_count);
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	auto before = high_resolution_clock::now();

	FileListLoader file_list(config.file_list);

    const int THREADS = 64;

	#pragma omp parallel num_threads(THREADS) 
	for (int i = omp_get_thread_num(); i < file_list.files.size(); i += omp_get_num_threads())
	{
        auto file_list_element = file_list.files[i];
		HashSortedArray hash_array;
		cout << "loading " << file_list_element.filename + config.allowed_kmers_postfix << endl;
		int kmer_len = DBSIO::load_dbs(file_list_element.filename + config.allowed_kmers_postfix, hash_array);
		auto centroid = build_centroid(file_list_element.filename, kmer_len, [&](hash_t kmer)
			{ 
				return has_hash(kmer, hash_array);
			});

		cout << centroid.seq_count << " " << centroid.kmers.size() << " " << file_list_element.filename + config.postfix << endl;
		DBSIO::save_centroid(file_list_element.filename + config.postfix, centroid, kmer_len);
	}

	cerr << "total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count() << endl;
}

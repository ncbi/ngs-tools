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
#include <list>
#include <stdexcept>
#include <iostream>
#include <chrono>
#include "omp_adapter.h"
#include "reader.h"
#include <unordered_map>

typedef uint64_t hash_t;

#include "log.h"
#include "dbs.h"
#include "fasta.h"
#include "hash.h"
#include "config_dump_kmers.h"
#include "seq_transform.h"
#include "acc_list_loader.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.12";

typedef std::unordered_map<hash_t, int> Kmers;

template <class Lambda>
void do_for_run(const string &filename, Lambda &&lambda)
{
    Reader::Params params;
    params.split_non_atgc = true;
    params.unaligned_only = false; //unaligned_only;
    auto reader = Reader::create(filename, params);

    std::vector<Reader::Fragment> chunk;
    bool done = false;
    while (!done) 
    {
        done = !reader->read_many(chunk);
        for (size_t seq_id = 0; seq_id < chunk.size(); ++seq_id) 
        {
//            auto& spotid = chunk[seq_id].spotid;
            auto& bases = chunk[seq_id].bases;
            lambda(bases);
        }
    }
}

Kmers get_kmer_set(const string &filename, int kmer_len)
{
    Kmers kmers;
    kmers.reserve(10000000);
    do_for_run(filename, [&](const string &seq)
    {
        Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
        {
			hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
            kmers[hash]++;
            return true;
        });
    });

    return kmers;
}

void print(const string &filename, const Kmers &kmers, int min_coverage, int kmer_len)
{
    ofstream f(filename);

    for (auto &k : kmers)
        if (k.second >= min_coverage)
            f << Hash<hash_t>::str_from_hash(k.first, kmer_len) << endl;
}

bool file_exists(const string &filename)
{
    ifstream f(filename);
    return f.good();
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	auto before = high_resolution_clock::now();

	AccListLoader acc_list(config.acc_list);

    const int THREADS = 48;
	#pragma omp parallel num_threads(THREADS) 
	for (int i = omp_get_thread_num(); i < acc_list.files.size(); i += omp_get_num_threads())
	{
        string out_file = acc_list.files[i] + ".kmers";
        if (file_exists(out_file))
            continue;

        #pragma omp critical
        {
            std::cerr << acc_list.files[i] << std::endl;
        }

		auto kmers = get_kmer_set(acc_list.files[i], config.kmer_len);
        print(out_file, kmers, config.min_coverage, config.kmer_len);
	}

	cerr << "total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count() << endl;
}

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

typedef uint64_t hash_t;

#include "log.h"
#include "dbs.h"
//#include "aligns_to_dbs_job.h"
#include "fasta.h"
#include "hash.h"
#include "config_fasta_contamination_multi.h"
#include "seq_transform.h"
#include "file_list_loader.h"

struct KmerTax : public DBS::KmerTaxMulti
{
    KmerTax() = default;
	KmerTax(hash_t kmer, const std::vector<int> &tax_ids) : DBS::KmerTaxMulti(kmer, tax_ids) { } // todo: remove constructor from KmerTax for faster loading ?

	bool operator < (const KmerTax &x) const // for binary search by hash
	{
		return kmer < x.kmer;
	}
};

using namespace std;
using namespace std::chrono;

const string VERSION = "0.10";

typedef int tax_t;
typedef vector<KmerTax> HashSortedArray;

const std::vector<int> &find_hash(hash_t hash, const std::vector<int> &default_value, const HashSortedArray &hash_array)
{
    auto first = hash_array.begin();
    auto last = hash_array.end();
    first = std::lower_bound(first, last, KmerTax(hash, default_value));
    return ((first == last) || (hash < first->kmer) ) ? default_value : first->tax_ids;
}

struct Hit
{
    int pos;
    tax_t tax_id;
    Hit(int pos, tax_t tax_id) : pos(pos), tax_id(tax_id){}
};

typedef std::list<Hit> Hits;

struct SeqContamination
{
    string desc;
    Hits hits;
    SeqContamination(const string &desc, const Hits &hits) : desc(desc), hits(hits){}
};

typedef std::list<SeqContamination> SeqContaminations;

struct Contamination
{
    string filename;
    SeqContaminations seqs;
    int total_hits;

    Contamination(const string &filename) : filename(filename), total_hits(0){}

    bool operator <(const Contamination &x) const
    {
        return total_hits > x.total_hits;
    }

    void add(const SeqContamination &seq)
    {
        seqs.push_back(seq);
        total_hits += seq.hits.size();
    }
};

Contamination check_for_contamination(const string &filename, int kmer_len, HashSortedArray &hash_array)
{
    Fasta fasta(filename);
    
    Contamination result(filename);

    std::vector<int> not_found;

    string seq;
    while (fasta.get_next_sequence(seq))
    {
        string desc = fasta.sequence_description();
        Hits hits;

        int pos = 0;
        Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
        {
			hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
            auto tax_ids = find_hash(hash, not_found, hash_array);
            for (auto tax_id : tax_ids)
                hits.push_back(Hit(pos, tax_id));

            pos++;
            return true;
        });

        if (!hits.empty())
            result.add(SeqContamination(desc, hits));
    }    

    return result;
}

void print(const Contamination &r)
{
    cout << r.filename << endl;
    cout << r.seqs.size() << endl;
    for (auto &s : r.seqs)
    {
        cout << s.desc << endl;
        cout << s.hits.size() << endl;
        for (auto &h : s.hits)
            cout << h.pos << '\t' << h.tax_id << endl;
    }
}

#define SORTED 0

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	auto before = high_resolution_clock::now();

	FileListLoader file_list(config.file_list);

	HashSortedArray hash_array;
	int kmer_len = DBSIO::load_dbsm(config.dbsm, hash_array);

    const int THREADS = 48;
#if SORTED
    list<Contamination> contaminations; // todo: list ?
#endif

	#pragma omp parallel num_threads(THREADS) 
	for (int i = omp_get_thread_num(); i < file_list.files.size(); i += omp_get_num_threads())
	{
        auto file_list_element = file_list.files[i];
		auto contamination = check_for_contamination(file_list_element.filename, kmer_len, hash_array);

        #pragma omp critical (read)
        {
		    cerr << contamination.total_hits << "\thits\t" << file_list_element.filename << endl; // must be stderr - stdout for output only
            if (contamination.total_hits > 0)
            {
#if SORTED
                contaminations.push_back(contamination);
#else
                print(contamination);
#endif
            }
        }
	}

#if SORTED
    contaminations.sort();

    for (auto &r : contaminations)
        print(r);
#endif

	cerr << "total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count() << endl;
}

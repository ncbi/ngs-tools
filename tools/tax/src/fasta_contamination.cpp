#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <stdexcept>
#include <iostream>
#include <chrono>
#include <omp.h>

typedef uint64_t hash_t;

#include "dbs.h"
//#include "aligns_to_dbs_job.h"
#include "fasta.h"
#include "hash.h"
#include "config_fasta_contamination.h"
#include "seq_transform.h"
#include "file_list_loader.h"

struct KmerTax : public DBS::KmerTax
{
	KmerTax(hash_t kmer = 0, int tax_id = 0) : DBS::KmerTax(kmer, tax_id){} // todo: remove constructor from KmerTax for faster loading ?

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

tax_t find_hash(hash_t hash, tax_t default_value, HashSortedArray &hash_array)
{
    auto first = hash_array.begin();
    auto last = hash_array.end();
    first = std::lower_bound(first, last, KmerTax(hash, 0));
    return ((first == last) || (hash < first->kmer) ) ? default_value : first->tax_id;
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

    string seq;
    while (fasta.get_next_sequence(seq))
    {
        string desc = fasta.sequence_description();
        Hits hits;

        int pos = 0;
        Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
        {
			hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
            auto tax_id = find_hash(hash, 0, hash_array);
            if (tax_id)
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

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	auto before = high_resolution_clock::now();

	FileListLoader file_list(config.file_list);

	HashSortedArray hash_array;
	int kmer_len = DBSIO::load_dbs(config.dbs, hash_array);

    const int THREADS = 48;
    list<Contamination> contaminations; // todo: list ?

	#pragma omp parallel num_threads(THREADS) 
	for (int i = omp_get_thread_num(); i < file_list.files.size(); i+=THREADS)
//	#pragma omp parallel for num_threads(THREADS) 
//    for (int i = 0; i < file_list.files.size(); i++)
	{
        auto file_list_element = file_list.files[i];
		auto contamination = check_for_contamination(file_list_element.filename, kmer_len, hash_array);

        #pragma omp critical (read)
        {
		    cerr << contamination.total_hits << "\thits\t" << file_list_element.filename << endl;
            if (contamination.total_hits > 0)
                contaminations.push_back(contamination);
//            if (contamination.total_hits > 0)
//                print(contamination);
        }
	}

//    std::sort(contaminations.begin(), contaminations.end());
//    cout << "sorted:" << endl;
    contaminations.sort();

    for (auto &r : contaminations)
        print(r);

	cerr << "total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count() << endl;
}

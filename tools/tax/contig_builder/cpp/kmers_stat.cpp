#include "config_kmers_stat.h"
#include <iostream>
#include <chrono>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include "fasta_reader.h"
#include "hash.h"
#include "log.h"
#include "seq_transform.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.1";

typedef uint32_t KmerCounter;

struct IndexBuilder// todo: multithreaded
{ 
	static const int KMER_LEN = 12;
	typedef uint32_t hash_t;

	//typedef vector<hash_t> Kmers;
	typedef std::unordered_map<hash_t, KmerCounter> KmerMap;
	KmerMap kmers;

	struct HashCount
	{
		hash_t hash;
		KmerCounter count;
		HashCount(hash_t hash, KmerCounter count) : hash(hash), count(count){}
		bool operator < (const HashCount &b) const { return count > b.count; }
	};


	vector<HashCount> kmers_stat;

	IndexBuilder()
	{
	}

	void add_reference(const string &seq)
	{
		Hash<hash_t>::for_all_hashes_do(seq, KMER_LEN, [&](hash_t hash)
			{
//				hash = seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN);
				kmers[hash]++;
				return true;
			});
	}

	void sort_kmers()
	{
		kmers_stat.reserve(kmers.size());
		for (auto const &rec : kmers)
			kmers_stat.push_back(HashCount(rec.first, rec.second));

		sort(kmers_stat.begin(), kmers_stat.end());
	}
};

void add_references(IndexBuilder &builder, FastaReader &fasta)
{
    Reader::Fragment fragment;
	while (fasta.read(&fragment))
	{
		if (fragment.bases.empty())
			return;
		builder.add_reference(fragment.bases);
	}
}

int main(int argc, char const *argv[])
{
    try
    {
		Config config(argc, argv);
		cerr << "kmers stat version " << VERSION << endl;

		print_current_time();
		auto before = high_resolution_clock::now();

		IndexBuilder builder;
		FastaReader fasta(config.reference);
		builder.kmers.reserve(fasta.file_size());
		cerr << "kmers reserved at (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;
		add_references(builder, fasta);
//		cerr << endl;
		cerr << "kmers found at (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;
		cerr << "kmers: " << builder.kmers.size() << endl;
		cout << "kmers: " << builder.kmers.size() << endl;
		builder.sort_kmers();
		cerr << "kmers sorted at (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;
//		builder.create_index();
//		cerr << "index built at (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;

//		save_index(config.out_file, builder.kmers, builder.buckets);
//		cerr << "total time (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;
//		int count = 0;
		for (auto &stat : builder.kmers_stat)
		{
			cout << stat.count << " " << Hash<IndexBuilder::hash_t>::str_from_hash(stat.hash, IndexBuilder::KMER_LEN) << endl;
			if (stat.count == 1)
				break;

//			count ++;
//			if (count > 100)
//				break;
		}

		exit(0); // dont want to wait for destructors
        return 0;
    }
    catch ( exception & x )
    {
        cerr << x.what() << endl;
//		cerr << "exit 3" << endl;
		return 3;
    }
    catch ( string & x )
    {
        cerr << x << endl;
//		cerr << "exit 4" << endl;
		return 4;
    }
    catch ( const char * x )
    {
        cerr << x << endl;
//		cerr << "exit 5" << endl;
		return 5;
    }
    catch ( ... )
    {
        cerr << "unknown exception" << endl;
//		cerr << "exit 6" << endl;
		return 6;
    }
}

#include "config_align_to_build_index_db.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include "fasta_reader.h"
#include "hash.h"
#include "log.h"
#include "seq_transform.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.1";

struct IndexBuilder// todo: multithreaded
{ 
	const int KMER_LEN = 32;
	typedef uint64_t hash_t;

	typedef vector<hash_t> Kmers;
	Kmers kmers;

	struct Bucket
	{
		unsigned int from, count;
		Bucket() : from(0), count(0){}
	};

	typedef vector<Bucket> Buckets;
	size_t BUCKETS_COUNT = size_t(1) << 32;

	Buckets buckets;

	IndexBuilder() : buckets(BUCKETS_COUNT)
	{
	}

	size_t bucket_of(hash_t hash)
	{
		return hash % BUCKETS_COUNT; // todo: optimize ?
	}

	void add_reference(const string &seq)
	{
		Hash<hash_t>::for_all_hashes_do(seq, KMER_LEN, [&](hash_t hash)
			{
				hash = seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN);
				kmers.push_back(hash);
				return true;
			});
	}

	void sort_kmers()
	{
		sort(kmers.begin(), kmers.end());
		cerr << "sorted!" << endl;
		kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end()); // todo: measure perf. probably can be removed
	}

	void create_index()
	{
		size_t i = 0;
		
		for (auto it = kmers.begin(); it != kmers.end(); it++, i++)
		{
			auto bucket = bucket_of(*it);
			if (!buckets[bucket].count)
				buckets[bucket].from = i;
			buckets[bucket].count++;
		}
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

void save_index(const string &filename, IndexBuilder::Kmers &kmers, IndexBuilder::Buckets &buckets)
{
	ofstream f(filename, ios::out | ios::binary);
	{
		size_t size = kmers.size();
		f.write((char*)&size, sizeof(size));
		f.write((char*)&kmers[0], kmers.size()*sizeof(kmers[0]));
	}

	{
		size_t size = buckets.size();
		f.write((char*)&size, sizeof(size));
		f.write((char*)&buckets[0], buckets.size()*sizeof(buckets[0]));
	}
}


int main(int argc, char const *argv[])
{
    try
    {
		Config config(argc, argv);
		cerr << "aligns_to_build_index_db version " << VERSION << endl;

		print_current_time();
		auto before = high_resolution_clock::now();

		IndexBuilder builder;
		FastaReader fasta(config.reference);
		builder.kmers.reserve(fasta.file_size());
		cerr << "kmers reserved at (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;
		add_references(builder, fasta);
//		cerr << endl;
		cerr << "kmers: " << builder.kmers.size() << endl;
		cerr << "kmers found at (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;
		builder.sort_kmers();
		cerr << "kmers sorted at (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;
		cerr << "kmers: " << builder.kmers.size() << endl;
		builder.create_index();
		cerr << "index built at (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;

		save_index(config.out_file, builder.kmers, builder.buckets);
		cerr << "total time (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;

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

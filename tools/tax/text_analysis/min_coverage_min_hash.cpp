#include <iostream>
#include <chrono>
#include <thread>
#include <list>
#include <iomanip>
#include <algorithm>
#include <map>
#include "config_min_coverage.h"
#include "text_loader_mt.h"
#include <omp.h>
#include <math.h>
#include <set>
#include "p_string.h"
#include <functional>

using namespace std;
using namespace std::chrono;

const string VERSION = "0.1";
#include "time.h"

int min_hash_pos(const string &s, int from, int kmer_len, int lookup_len)
{
	std::hash<std::string> hash_fn;
	size_t min_hash = hash_fn(string(&s[from], kmer_len));
	int best_pos = from;
	for (int pos = from + 1; pos < from + lookup_len; pos++)
	{
		auto h = hash_fn(string(&s[pos], kmer_len));
		if (h < min_hash)
		{
			min_hash = h;
			best_pos = pos;
		}
	}

	return best_pos;
}

template <class OnKmerFound, class OnSeqCompressed>
void compress(TextLoaderSTNoStore &seq_loader, int kmer_len, int lookup_len, OnKmerFound &&kmer_found, OnSeqCompressed &&seq_compressed)
{
	string seq;
	while (seq_loader.load_next_sequence(seq))
	{
		for (int pos = 0; pos < int(seq.length()) - kmer_len - lookup_len; pos += kmer_len + lookup_len)
			kmer_found(seq, min_hash_pos(seq, pos, kmer_len, lookup_len));

		seq_compressed(seq);
	}
}

int main(int argc, char const *argv[])
{
    try
    {
		ConfigMinCoverage config(argc, argv);
		cerr << "min_coverage_min_hash version " << VERSION << endl;

		print_current_time();
		auto before = high_resolution_clock::now();

		TextLoaderSTNoStore seq_loader(config.reference);
		std::set<string> kmers;

		compress(seq_loader, config.kmer_len, config.lookup_len, 
			[&](const string &seq, size_t pos)
			{
//				print_string(&seq[pos], config.kmer_len);
//				cout << endl;
				kmers.insert(string(&seq[pos], config.kmer_len));
			},
			[](const string &seq)
			{
				cerr << ".";
			}
		);

//		cout << "-------" << endl;
		for (auto &s : kmers)
			cout << s << endl;

		cerr << "total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;

		exit(0); // dont want to wait for KmerMaps destructors
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

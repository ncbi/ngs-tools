#include <iostream>
#include <chrono>
#include <thread>
#include <list>
#include <algorithm>
#include <map>
#include "config_kmer_distribution.h"
#include "text_loader_mt.h"
#include <omp.h>

using namespace std;
using namespace std::chrono;

const string VERSION = "0.1";

int global_kmer_len = 0; // todo: think about better solution. 

int string_compare(const char *a, const char *b, unsigned int len)
{
	for (unsigned int i=0; i<len; i++)
	{
		int diff = int(a[i]) - int(b[i]);
		if (diff != 0)
			return diff;
	}

	return 0;
}

void print_string(const char *s, int len)
{
	for (int i=0; i<len; i++)
		cout << s[i];
}


void print_current_time()
{
	auto t = std::time(nullptr);
	auto timeinfo = std::localtime(&t);
	const int BUFFER_SIZE = 256;
	char buffer[BUFFER_SIZE];
	strftime(buffer, BUFFER_SIZE, "%m/%d/%Y %H:%M:%S", timeinfo);
	cerr << "time is " << buffer << endl;
}


struct KmerDict
{
	typedef unsigned int Count;

	struct Key
	{
		const char *s;
		Key(const char *s) : s(s){}
		bool operator < (const Key &b) const
		{
			return string_compare(s, b.s, global_kmer_len) < 0;
		}

		//bool operator == (const Key &b) const
		//{
		//	return string_compare(s, b.s, global_kmer_len) == 0;
		//}
	};

	std::map<Key, Count> storage;

	void add_kmer(const char *s, Count amount = 1)
	{
		//cerr << "adding string"; // << (void*)s;
		//print_string(s, global_kmer_len);
		//cerr << endl;

		auto key = Key(s);
		auto &c = storage[key];
		c += amount;
	}
};

void build_kmers(KmerDict &kmer_dict, TextLoaderMT &text_loader)
{
	while (true)
	{
		string *seq = text_loader.load_next_sequence();
		if (seq == nullptr)
			return;

		for (int i = 0; i < int(seq->size()) - global_kmer_len + 1; i++)
			kmer_dict.add_kmer(&((*seq)[i]));
	}
}

void merge(KmerDict &target_dict, const KmerDict &kmer_dict, unsigned int min_coverage)
{
	for (auto &map_element: kmer_dict.storage)
		if (map_element.second > min_coverage) // it's just optimization - not completely correct
			target_dict.add_kmer(map_element.first.s, map_element.second);
}

void print_result_kmer(KmerDict::Key key, KmerDict::Count count)
{
	print_string(key.s, global_kmer_len);
	cout << "|" << count << endl;
}

void print_result(KmerDict &dict, unsigned int min_coverage)
{
	struct KeyCount
	{
		KmerDict::Key key;
		KmerDict::Count count;
		KeyCount(KmerDict::Key key, KmerDict::Count count) : key(key), count(count){}

		bool operator < (const KeyCount &x) const
		{
			return count > x.count;
		}
	};

	vector<KeyCount> v;
	v.reserve(dict.storage.size());
	for (auto &map_element: dict.storage)
		if (map_element.second >= min_coverage)
		{
//			cout << (void*)map_element.first.s;
			v.push_back(KeyCount(map_element.first, map_element.second));
		}

	std::sort(v.begin(), v.end());

	for (auto &x : v)
		print_result_kmer(x.key, x.count);
}

const int THREADS = 16;

int main(int argc, char const *argv[])
{
    try
    {
		ConfigKmerDistribution config(argc, argv);
		cerr << "kmer_distibution version " << VERSION << endl;

		print_current_time();
		auto before = high_resolution_clock::now();

		TextLoaderMT seq_loader(config.reference);
		vector<KmerDict> kmer_dicts(THREADS);
//		for (int i=0; i<THREADS; i++)
//			kmer_dicts.push_back(KmerDict(config.kmer_len));
		global_kmer_len = config.kmer_len; // todo: think

		#pragma omp parallel num_threads(THREADS)
		{
			auto i = omp_get_thread_num();
			build_kmers(kmer_dicts[i], seq_loader);
		}

		cerr << "merging" << endl;
		cerr << "spent time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;
		auto &final_dict = kmer_dicts[0];
		for (int i=1; i<THREADS; i++)
		{
			merge(final_dict, kmer_dicts[i], config.min_coverage);
			cerr << ".";
		}
		cerr << endl;

		cerr << "printing results" << endl;
		print_result(final_dict, config.min_coverage);

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

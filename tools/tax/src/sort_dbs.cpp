#include "config_sort_dbs.h"
#include <iostream>
#include <vector>
#include <list>
#include <fstream>
#include <algorithm>

using namespace std;

const string VERSION = "0.12";

typedef uint64_t hash_t;

#include "dbs.h"

typedef DBS::Kmers Kmers;
typedef vector<hash_t> Hashes;

bool split_by_tax_less(const DBS::KmerTax &a, const DBS::KmerTax &b)
{
	if (a.tax_id == b.tax_id)
		return a.kmer < b.kmer;

	return a.tax_id < b.tax_id;
}

void to_hashes(const Kmers &kmers, Hashes &hashes)
{
	hashes.clear();
	hashes.resize(kmers.size());
	for (size_t i = 0; i < kmers.size(); i++)
		hashes[i] = kmers[i].kmer;
};

struct Annot
{
	int tax_id;
	size_t count;

	Annot(int tax_id, size_t count) : tax_id(tax_id), count(count){}
};

typedef list<Annot> Annotation;

void save_annotation(const string &filename, const Annotation &annotation)
{
	ofstream f(filename); //, std::ios::text| std::ios::out);
	f.flush();
	for (auto &a : annotation)
		f << a.tax_id << '\t' << a.count << endl;
}

Annotation get_annotation(const Kmers &kmers)
{
	Annotation a;
	if (kmers.empty())
		return a;

	size_t start = 0;
	auto current_tax = kmers[0].tax_id;
	for (size_t i = 1; i < kmers.size(); i++)
		if (kmers[i].tax_id != current_tax)
		{
			a.push_back(Annot(current_tax, i - start));
			start = i;
			current_tax = kmers[i].tax_id;
			if (current_tax <= 0)
				throw "invalid taxonomy";
		}

	a.push_back(Annot(current_tax, kmers.size() - start));
	return a;
}

int main(int argc, char const *argv[])
{
    try
    {
		Config config(argc, argv);
		cerr << "sort_dbs version " << VERSION << endl;

		Kmers kmers;
		auto kmer_len = DBS::load_dbs(config.input_filename, kmers);
		std::sort(kmers.begin(), kmers.end(), split_by_tax_less);
		Hashes hashes;
		to_hashes(kmers, hashes);
		DBS::save_dbs(config.out_filename, hashes, kmer_len);
		save_annotation(config.out_filename + ".annotation", get_annotation(kmers));

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

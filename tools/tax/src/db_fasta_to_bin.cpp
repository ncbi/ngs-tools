#include "config_db_fasta_to_bin.h"
#include "text_loader_mt.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include "hash.h"

using namespace std;
//using namespace std::chrono;

//const int KMER_LEN = 32;
typedef uint64_t hash_t;

#include "dbs.h"

const string VERSION = "0.22";

char complement(char ch)
{
	switch (ch)
	{
		case 'A': return 'T';
		case 'T': return 'A';
		case 'C': return 'G';
		case 'G': return 'C';
	};

	throw string("bad letter to reverse: ") + string(1, ch);
}

string reverse_complement(string s)
{
	std::reverse(s.begin(), s.end());
	for (char &ch : s)
		ch = complement(ch);

	return s;
}

hash_t hash_of(const string &s)
{
	if (!s.size())
		throw string("invalid string len: 0 ");

	return Hash<hash_t>::hash_of(s);
}

void process_without_taxonomy(const string &fasta_db, const string &out_file)
{
	cout << "process without taxonomy info" << endl;
	TextLoaderSTNoStore loader(fasta_db);
	vector<hash_t> kmers;

	int kmer_len = 0;
	{
		string seq;
		while (loader.load_next_sequence(seq))
		{
			kmers.push_back(std::min( hash_of(seq), hash_of(reverse_complement(seq))) );
			if (!kmer_len)
				kmer_len = seq.length();

			if (seq.length() != kmer_len)
				throw "seq.length() != kmer_len";
		}
	}

	sort(kmers.begin(), kmers.end());
	DBS::save_dbs(out_file, kmers, kmer_len);
}

typedef DBS::KmerTax KmerTax;

bool kmer_less(const KmerTax &a, const KmerTax &b)
{
	return a.kmer < b.kmer;
}

void process_with_taxonomy(const string &fasta_db, const string &out_file)
{
	cout << "process with taxonomy info" << endl;

	FastaWithTaxonomyLoader loader(fasta_db);
	vector<KmerTax> kmers;
	int kmer_len = 0;

	{
		string seq;
		int tax_id;
		while (loader.load_next_sequence(seq, tax_id))
		{
			kmers.push_back(KmerTax(std::min( hash_of(seq), hash_of(reverse_complement(seq))), tax_id) );
			if (!kmer_len)
				kmer_len = seq.length();

			if (seq.length() != kmer_len)
				throw "seq.length() != kmer_len";
		}
	}

	sort(kmers.begin(), kmers.end(), kmer_less);
	DBS::save_dbs(out_file, kmers, kmer_len);
}

bool has_taxonomy_info(const string &filename)
{
	ifstream f(filename);
	if (f.fail() || f.eof())
		throw string("cannot open file ") + filename;

	string seq1, seq2;
	f >> seq1;
	f >> seq2;

	return seq1.length() != seq2.length();
}

int main(int argc, char const *argv[])
{
    try
    {
		Config config(argc, argv);
		cerr << "db_fasta_to_bin version " << VERSION << endl;

		if (has_taxonomy_info(config.fasta_db))
			process_with_taxonomy(config.fasta_db, config.out_file);
		else
			process_without_taxonomy(config.fasta_db, config.out_file);

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

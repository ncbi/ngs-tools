#include <iostream>
#include <fstream>

#include "config_filter_db.h"
using namespace std;

const string VERSION = "0.11";

void dont_pass(const string &kmer, unsigned int tax_id)
{
	cerr << kmer << endl;
}

void pass(const string &kmer, unsigned int tax_id)
{
	cout << kmer << '\t' << tax_id << endl;
}

void pass(const string &kmer)
{
	cout << kmer << endl;
}

int predicted(const string &kmer)
{
	int pred = 0;
	for (int i = 4; i<int(kmer.length()); i++)
		if (kmer[i] == kmer[i - 1] && kmer[i - 1] == kmer[i - 2] && kmer[i - 2] == kmer[i - 3] && kmer[i - 3] == kmer[i - 4])
			pred++;
		else if (i >= 8 && kmer[i] == kmer[i - 2] && kmer[i - 2] == kmer[i - 4] && kmer[i - 4] == kmer[i - 6] && kmer[i - 6] == kmer[i - 8])
			pred++;

	return pred;
}

bool bad_tax(unsigned int tax_id, unsigned int config_only_tax_id)
{
	const int CELLULAR_ORGANISMS = 131567;
	if (tax_id == CELLULAR_ORGANISMS)
		return true;

	return config_only_tax_id && (tax_id != config_only_tax_id);
}

int main(int argc, char const *argv[])
{
    try
    {
		ConfigFilterDB config(argc, argv);
		cerr << "filter_db version " << VERSION << endl;
		ifstream f(config.input_file);
		if (f.fail())
			throw "cannot open input file";

		if (config.only_tax)
			cerr << "keep only tax " << config.only_tax << endl;

		string kmer;
		unsigned int tax_id;
		while (!f.eof())
		{
			f >> kmer;
			if (!kmer.length() || f.eof())
				break;

			f >> tax_id;
			if (bad_tax(tax_id, config.only_tax))
			{
				dont_pass(kmer, tax_id);
				continue;
			}

			auto score = predicted(kmer);
			const int MIN_SCORE = 13 * kmer.length()/32; // const 13 was designed for 32 bp kmers
			if (score >= MIN_SCORE)
				dont_pass(kmer, tax_id);
			else
				if (config.only_tax)
					pass(kmer);
				else
					pass(kmer, tax_id);
		}

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

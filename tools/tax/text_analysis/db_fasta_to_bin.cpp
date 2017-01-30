#include "config_db_fasta_to_bin.h"
#include "text_loader_mt.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include "hash.h"

using namespace std;
//using namespace std::chrono;

const string VERSION = "0.1";

const int KMER_LEN = 32;
typedef uint64_t hash_t;

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
	if (s.size() != KMER_LEN)
		throw string("invalid string len: ") + s;

	return Hash<hash_t>::hash_of(s);
}

int main(int argc, char const *argv[])
{
    try
    {
		Config config(argc, argv);
		cerr << "db_fasta_to_bin version " << VERSION << endl;

		TextLoaderSTNoStore loader(config.fasta_db);
		vector<hash_t> kmers;

		{
			string seq;
			while (loader.load_next_sequence(seq))
				kmers.push_back(std::min( hash_of(seq), hash_of(reverse_complement(seq))) );
		}

		sort(kmers.begin(), kmers.end());
		cout << "writing " << kmers.size() << " kmers" << endl;

		{
			ofstream f(config.out_file);
			f.flush();
			size_t kmers_size = kmers.size();
			f.write((char*)&kmers_size, sizeof(kmers_size));
			f.write((char*)&kmers[0], kmers.size()*sizeof(kmers[0]));
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

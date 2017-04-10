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

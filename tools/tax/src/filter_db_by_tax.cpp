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
#include <stdexcept>

#include "log.h"
#include "config_filter_db_by_tax.h"

using namespace std;

const string VERSION = "0.10";

void pass(const string &kmer)
{
	cout << kmer << endl;
}

bool bad_tax(unsigned int tax_id, unsigned int config_only_tax_id)
{
	return config_only_tax_id && (tax_id != config_only_tax_id);
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	LOG("filter_db_by_tax version " << VERSION);
		
	ifstream f(config.input_file);
	if (f.fail())
		throw std::runtime_error("cannot open input file");

	if (config.only_tax)
		LOG("keep only tax " << config.only_tax);

	string kmer;
	unsigned int tax_id;

	while (!f.eof())
	{
		f >> kmer;
		if (!kmer.length() || f.eof())
			break;

		f >> tax_id;
		if (bad_tax(tax_id, config.only_tax))
			continue;

		pass(kmer);
	}

    return 0;
}

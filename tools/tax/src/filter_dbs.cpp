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
#include <chrono>
#include <thread>
#include <array>
#include "kmers.h"
#include "filter_db.h"
#include "dbs.h"
#include "hash.h"
#include "config_filter_dbs.h"

using namespace std;

bool bad_id(int tax_id)
{
	const int ROOT = 1;
	return tax_id == 0 || tax_id == ROOT;
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	std::vector<DBS::KmerTax> kmers;
	auto kmer_len = DBSIO::load_dbs(config.in_file, kmers);

	size_t good = 0, bad_ids = 0, low_comp = 0;
	for (size_t i = 0; i < kmers.size(); i++)
		if (bad_id(kmers[i].tax_id))
			bad_ids++;
		else if (FilterDB::low_complexity(Hash<hash_t>::str_from_hash(kmers[i].kmer, kmer_len)))
			low_comp++;
		else
			good++;

	cout << "out of " << kmers.size() << endl;
	cout << "good: " << good << endl;
	cout << "bad id: " << bad_ids << endl;
	cout << "low complexity : " << low_comp << endl;

	std::ofstream f(config.out_file);
	DBSIO::DBSHeader header(kmer_len);
	IO::write(f, header);

	IO::write(f, good);
	size_t written = 0;
	for (size_t i = 0; i < kmers.size(); i++)
		if (!bad_id(kmers[i].tax_id) && !FilterDB::low_complexity(Hash<hash_t>::str_from_hash(kmers[i].kmer, kmer_len)))
		{
			IO::write(f, kmers[i]);
			written++;
		}

	if (written != good)
		throw std::runtime_error("written != good");
}


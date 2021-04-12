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
#include "config_db_tax_id_to_dbs.h"

using namespace std;

bool bad(hash_t kmer, int tax_id)
{
	return tax_id == 0; // todo: complexity check?
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	std::vector<hash_t> kmers;
	auto kmer_len = DBSIO::load_dbs(config.db_file, kmers);

	std::vector<int> tax_ids;

	{
		std::ifstream f(config.tax_id_file);
		IO::load_vector(f, tax_ids);
		if (kmers.size() != tax_ids.size())
			throw std::runtime_error("kmers.size() != tax_ids.size()");
	}

	size_t good = 0;
	for (size_t i = 0; i < kmers.size(); i++)
		if (bad(kmers[i], tax_ids[i]))
			tax_ids[i] = 0;
		else
			good++;

	std::ofstream f(config.out_dbs_file);
	DBSIO::DBSHeader header(kmer_len);
	IO::write(f, header);

	IO::write(f, good);
	size_t written = 0;
	for (size_t i = 0; i < kmers.size(); i++)
		if (tax_ids[i] != 0)
		{
			IO::write(f, DBS::KmerTax(kmers[i], tax_ids[i]));
			written++;
		}

	if (written != good)
		throw std::runtime_error("written != good");
}


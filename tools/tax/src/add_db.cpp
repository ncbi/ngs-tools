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
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include "dbs.h"
#include "hash.h"
#include "config_add_db.h"

bool in_db(hash_t hash, const std::vector<hash_t> &db, size_t lookup_len)
{
	return std::binary_search(db.begin(), db.begin() + lookup_len, hash);
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	std::vector<hash_t> db_a, db_b;
	int kmer_len_a = DBSIO::load_dbs(config.db_a, db_a);
	int kmer_len_b = DBSIO::load_dbs(config.db_b, db_b);

    if (kmer_len_a != kmer_len_b)
    {
        std::cerr << "kmer_len_a != kmer_len_b" << std::endl;
        return 1;
    }

    size_t original_a_size = db_a.size();
    db_a.reserve(db_a.size() + db_b.size());

	for (auto kmer : db_b)
        if (!in_db(kmer, db_a, original_a_size))
    		db_a.push_back(kmer);

	std::sort(db_a.begin(), db_a.end());
	DBSIO::save_dbs(config.db_c, db_a, kmer_len_a);
}

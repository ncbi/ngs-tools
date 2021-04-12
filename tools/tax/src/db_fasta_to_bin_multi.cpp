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

#include "config_db_fasta_to_bin_multi.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include "hash.h"
#include "seq_transform.h"
#include "log.h"
#include <sstream>

using namespace std;
typedef uint64_t hash_t;

#include "dbs.h"

const string VERSION = "0.10";

void fail(const std::string &message)
{
	LOG(message);
	throw std::runtime_error(message);
}

hash_t hash_of(const string &s)
{
	if (!s.size())
		throw std::runtime_error("invalid string len: 0 ");

	return Hash<hash_t>::hash_of(s);
}

typedef DBS::KmerTaxMulti KmerTaxMulti;

bool kmer_less(const KmerTaxMulti &a, const KmerTaxMulti &b)
{
	return a.kmer < b.kmer;
}

void process_with_taxonomy(const string &filename, const string &out_file)
{
	vector<KmerTaxMulti> kmers;
	int kmer_len = 0;

	ifstream f(filename);
	if (f.fail() || f.eof())
		fail(string("cannot open file ") + filename);

	string kmer, line;
	while (!f.eof())
	{
        std::getline(f, line);
        if (line.empty())
            if (!f.eof())
                fail("empty line in the middle of kmers file");
            else
                break;

        istringstream iss(line);
		iss >> kmer;
		if (kmer.empty())
        {
            if (!f.eof())
                fail(string("error reading ") + filename);
			break;
        }

		if (!kmer_len)
			kmer_len = kmer.length(); // len of first kmer

		if (kmer.length() != kmer_len)
			throw std::runtime_error("seq.length() != kmer_len");

        auto hash = hash_of(kmer);
        if (seq_transform<hash_t>::min_hash_variant(hash, kmer_len) != hash)
			throw std::runtime_error("not normalized data : min hash variant != hash");

        if (iss.eof())
            fail("no tax id specified for kmer");

        vector<int> taxes;
        while (!iss.eof())
        {
    		int tax_id = -1;
		    iss >> tax_id;
		    if (iss.fail())
			    fail("bad tax id while loading kmers");

            taxes.push_back(tax_id);
        }

        kmers.push_back(KmerTaxMulti(hash, taxes));
	}

	sort(kmers.begin(), kmers.end(), kmer_less);
	DBSIO::save_dbsm(out_file, kmers, kmer_len);
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);
	LOG("db_fasta_to_bin_multi version " << VERSION);

	process_with_taxonomy(config.fasta_db, config.out_file);

    return 0;
}

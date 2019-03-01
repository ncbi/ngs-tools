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

#include "config_db_fasta_to_bin.h"
#include "text_loader.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include "hash.h"
#include "seq_transform.h"
#include "log.h"

using namespace std;
typedef uint64_t hash_t;

#include "dbs.h"

const string VERSION = "0.22";

string reverse_complement(string s) // yes, by value
{
    seq_transform_actg::to_rev_complement(s);
    return s;
}

hash_t hash_of(const string &s)
{
	if (!s.size())
		throw std::runtime_error("invalid string len: 0 ");

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
			kmers.push_back(std::min( hash_of(seq), hash_of(reverse_complement(seq))));
			if (!kmer_len)
				kmer_len = seq.length();

			if (seq.length() != kmer_len)
				throw std::runtime_error("seq.length() != kmer_len");
		}
	}

	sort(kmers.begin(), kmers.end());
	DBSIO::save_dbs(out_file, kmers, kmer_len);
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
			kmers.push_back(KmerTax(std::min( hash_of(seq), hash_of(reverse_complement(seq))), tax_id));
			if (!kmer_len)
				kmer_len = seq.length();

			if (seq.length() != kmer_len)
				throw std::runtime_error("seq.length() != kmer_len");
		}
	}

	sort(kmers.begin(), kmers.end(), kmer_less);
	DBSIO::save_dbs(out_file, kmers, kmer_len);
}

bool has_taxonomy_info(const string &filename)
{
	ifstream f(filename);
	if (f.fail() || f.eof())
		throw std::runtime_error(string("cannot open file ") + filename);

	string seq1, seq2;
	f >> seq1;
	f >> seq2;

	return seq1.length() != seq2.length();
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);
	LOG("db_fasta_to_bin version " << VERSION);

	if (has_taxonomy_info(config.fasta_db))
		process_with_taxonomy(config.fasta_db, config.out_file);
	else
		process_without_taxonomy(config.fasta_db, config.out_file);

    return 0;
}

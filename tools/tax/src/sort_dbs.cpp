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

#include "config_sort_dbs.h"
#include <iostream>
#include <vector>
#include <list>
#include <fstream>
#include <algorithm>
#include <stdint.h>

#include "log.h"

using namespace std;

const string VERSION = "0.12";

typedef uint64_t hash_t;

#include "dbs.h"

typedef DBS::Kmers Kmers;
typedef vector<hash_t> Hashes;

bool split_by_tax_less(const DBS::KmerTax &a, const DBS::KmerTax &b)
{
	if (a.tax_id == b.tax_id)
		return a.kmer < b.kmer;

	return a.tax_id < b.tax_id;
}

void to_hashes(const Kmers &kmers, Hashes &hashes)
{
	hashes.clear();
	hashes.resize(kmers.size());
	for (size_t i = 0; i < kmers.size(); i++)
		hashes[i] = kmers[i].kmer;
};

struct Annot
{
	int tax_id;
	size_t count;

	Annot(int tax_id, size_t count) : tax_id(tax_id), count(count){}
};

typedef list<Annot> Annotation;

void save_annotation(const string &filename, const Annotation &annotation)
{
	ofstream f(filename);
	for (auto &a : annotation)
		f << a.tax_id << '\t' << a.count << endl;
}

Annotation get_annotation(const Kmers &kmers)
{
	Annotation a;
	if (kmers.empty())
		return a;

	size_t start = 0;
	auto current_tax = kmers[0].tax_id;
	for (size_t i = 1; i < kmers.size(); i++)
		if (kmers[i].tax_id != current_tax)
		{
			a.push_back(Annot(current_tax, i - start));
			start = i;
			current_tax = kmers[i].tax_id;
			if (current_tax <= 0)
				throw std::runtime_error("invalid taxonomy");
		}

	a.push_back(Annot(current_tax, kmers.size() - start));
	return a;
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);
	LOG("sort_dbs version " << VERSION);

	Kmers kmers;
	auto kmer_len = DBSIO::load_dbs(config.input_filename, kmers);
	std::sort(kmers.begin(), kmers.end(), split_by_tax_less);
	Hashes hashes;
	to_hashes(kmers, hashes);
	DBSIO::save_dbs(config.out_filename, hashes, kmer_len);
	save_annotation(config.out_filename + ".annotation", get_annotation(kmers));

    return 0;
}

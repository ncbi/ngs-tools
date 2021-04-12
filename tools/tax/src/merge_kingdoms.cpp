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
#include <utility>
#include <map>
#include "kmers.h"
#include "dbs.h"
#include "tax_id_tree.h"
#include "config_merge_kingdoms.h"

using namespace std;

int merge_ids(int higher, int lower, int min_higher_level, const TaxIdTree &tax_id_tree)
{
	if (tax_id_tree.level(higher) > min_higher_level)
		return higher;

	return tax_id_tree.consensus_of(higher, lower);		
}

void print_stat(map< pair<int, int>, int> &stat, int min_higher_level, const TaxIdTree &tax_id_tree)
{
	for (auto &it : stat)
		cout << it.second << ": " << it.first.first << " " << it.first.second << " -> " << merge_ids(it.first.first, it.first.second, min_higher_level, tax_id_tree) << endl;
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	TaxIdTree tax_id_tree;
	TaxIdTreeLoader::load_tax_id_tree(tax_id_tree, config.tax_parents_file);

	std::vector<DBS::KmerTax> a, b;
	auto a_kmer_len = DBSIO::load_dbs(config.higher_file, a);
	auto b_kmer_len = DBSIO::load_dbs(config.lower_file, b);

	if (a_kmer_len != b_kmer_len)
		throw std::runtime_error("a_kmer_len != b_kmer_len");

	std::vector<DBS::KmerTax> merged;
	merged.reserve(a.size() + b.size());

	map< pair<int, int>, int> stat;

	size_t ai = 0, bi = 0;
	while (ai < a.size() && bi < b.size())
	{
		if (a[ai].kmer < b[bi].kmer)
		{
			merged.push_back(a[ai]);
			ai++;
		}
		else if (a[ai].kmer > b[bi].kmer)
		{
			merged.push_back(b[bi]);
			bi++;
		}
		else
			{
				auto tax_id = merge_ids(a[ai].tax_id, b[bi].tax_id, config.min_higher_level, tax_id_tree);
				merged.push_back(DBS::KmerTax(a[ai].kmer, tax_id));
				stat[make_pair(a[ai].tax_id, b[bi].tax_id)]++;
				ai++;
				bi++;
			}
	}

	for (; ai < a.size(); ai++)
		merged.push_back(a[ai]);
	for (; bi < b.size(); bi++)
		merged.push_back(b[bi]);

	cout << "     a size: " << a.size() << endl;
	cout << "     b size: " << b.size() << endl;
	cout << "merged size: " << merged.size() << endl;
	print_stat(stat, config.min_higher_level, tax_id_tree);

	DBSIO::save_dbs(config.out_file, merged, a_kmer_len);
}


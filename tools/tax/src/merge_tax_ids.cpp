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
#include <vector>
#include <chrono>
#include "config_merge_tax_ids.h"
#include "file_list_loader.h"
#include "io.h"
#include "tax_id_tree.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.10";

typedef vector<int> TaxIds;

TaxIds load_tax_ids(const string &filename)
{
	TaxIds tax_ids;
	std::ifstream f(filename);
	IO::load_vector(f, tax_ids);
	return tax_ids;
}

void merge_ids(TaxIds &target, const TaxIds &b, const TaxIdTree &tax_id_tree)
{
	if (target.size() != b.size())
		throw std::runtime_error("tax ids have different sizes and cannot be merged");
	
	for (size_t i = 0; i < target.size(); i++)
		if (b[i] != 0)
		{
			auto &t = target[i];
			if (t == 0)
				t = b[i];
			else 
				t = tax_id_tree.consensus_of(t, b[i]);
		}
}

int main(int argc, char const *argv[])
{
	cout << "merge_tax_ids version " << VERSION << endl;
	Config config(argc, argv);

	auto before = high_resolution_clock::now();

	TaxIdTree tax_id_tree;
	TaxIdTreeLoader::load_tax_id_tree(tax_id_tree, config.tax_parents_file);

	FileListLoader file_list(config.file_list);

	TaxIds tax_ids;
	for (auto &file_list_element : file_list.files)
	{
		cout << file_list_element.filename << endl;
		auto loaded_ids = load_tax_ids(file_list_element.filename);
		if (tax_ids.empty())
			tax_ids = loaded_ids;
		else
			merge_ids(tax_ids, loaded_ids, tax_id_tree);
	}

	std::ofstream f(config.out_file);
	IO::save_vector(f, tax_ids);

	cout << "total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count() << endl;
}

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
#include "config_merge_db.h"
#include "file_list_loader.h"
#include "dbs.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.14";

template <class C>
struct HashTable
{
	std::vector<C> storage;
	size_t inserted = 0;
	HashTable(size_t storage_size) : storage(storage_size){}

	void insert(C x)
	{
		if (x == 0)
			return;

		for (size_t pos = x % storage.size(); storage[pos] != x ; pos = next(pos))
			if (storage[pos] == 0)
			{
				storage[pos] = x;
				inserted++;
				if (inserted > storage.size() * 2 / 3)
					throw std::runtime_error("HashTable inserted > storage_size * 2 / 3");
				break;
			}
	}

	size_t next(size_t pos) const
	{
		pos++;
		if (pos < storage.size())
			return pos;
		if (pos > storage.size())
			throw std::runtime_error("HashTable pos > storage.size()");

		return 0;		
	}

};

int main(int argc, char const *argv[])
{
	cout << "merge_db version " << VERSION << endl;
	Config config(argc, argv);

	auto before = high_resolution_clock::now();

	FileListLoader file_list(config.file_list);
	HashTable<hash_t> kmers(config.expected_size * 2);

	size_t kmer_len = 0;

	vector<hash_t> loaded_kmers;
	for (auto &file_list_element : file_list.files)
	{
		string in_file = file_list_element.filename + config.in_ext;
		cout << in_file << endl;
		auto loaded_kmer_len = DBSIO::load_dbs(in_file, loaded_kmers);
		if (kmer_len != 0 && loaded_kmer_len != kmer_len)
			throw std::runtime_error("loaded_kmer_len != kmer_len");

		if (loaded_kmers.empty())
		{
			cout << "skip empty" << endl;
			continue;
		}

		auto before_kmers = kmers.inserted;
		kmer_len = loaded_kmer_len;
		for (auto kmer : loaded_kmers)
			kmers.insert(kmer);

		auto inserted = kmers.inserted - before_kmers;
		cout << "inserted " << inserted << " of " << loaded_kmers.size() << " = " << int(100.0 * inserted / loaded_kmers.size()) << " %, fullness = " << int(kmers.inserted * 100.0 / kmers.storage.size()) << " %" << endl;
	}

	std::sort(kmers.storage.begin(), kmers.storage.end());
	cout << "saving to " << config.out_file << " kmers: " << kmers.inserted << endl;
	DBSIO::save_dbs(config.out_file, kmers.storage, kmer_len, kmers.storage.size() - kmers.inserted);

	cout << "total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count() << endl;
}

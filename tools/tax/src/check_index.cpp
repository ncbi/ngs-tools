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
#include "omp_adapter.h"
#include "check_index.h"
#include "kmers.h"
#include "kmer_io.h"
#include "ready_seq.h"
#include "filename_meta.h"
#include "file_list_loader.h"

#include "ready_seq.h"
#include <sstream>

#include "tax_id_tree.h"
#include "config_check_index.h"
#include "filename_meta.h"

using namespace std;
using namespace std::chrono;


void fail(const std::string &message)
{
	LOG(message);
	throw std::runtime_error(message);
}

int load_kmers(Kmers &kmers, const string &filename)
{
	int kmer_len = 0;

	ifstream f(filename);
	if (f.fail() || f.eof())
		fail(string("cannot open file ") + filename);

	string kmer, line;
	while (!f.eof())
	{
		f >> kmer;
		if (kmer.empty())
        {
            if (!f.eof())
               fail(string("error reading ") + filename);
			break;
        }                

		if (!kmer_len)
			kmer_len = kmer.length(); // len of first kmer

		if (kmer.length() != kmer_len)
			fail("wrong kmer len for " + kmer);

		int tax_id;
		f >> tax_id;
		if (!tax_id)
			fail("bad tax id while loading kmers");

		{
			hash_t hash = Hash<hash_t>::hash_of(kmer);
			hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
			kmers.add_kmer(hash, tax_id);
		}
	}

	return kmer_len;
}

int main(int argc, char const *argv[])
{
	ConfigCheckIndex config(argc, argv);
	LOG("check_index version 0.10 ");

	auto before = high_resolution_clock::now();

	FileListLoader file_list(config.file_list);

	TaxIdTree tax_id_tree;
	TaxIdTreeLoader::load_tax_id_tree(tax_id_tree, config.tax_parents_file);

	Kmers kmers(tax_id_tree);
	int kmer_len = load_kmers(kmers, config.kmers_file);
	LOG("kmer len: " << kmer_len);
	LOG(kmers.storage.size() << " kmers loaded");

	size_t total_size = 0;
	for (auto &file_list_element : file_list.files)
	{
		auto tax_id = FilenameMeta::tax_id_from(file_list_element.filename);
		LOG(file_list_element.filesize << "\t" << tax_id << "\t" << file_list_element.filename);
		total_size += CheckIndex<Kmers>::check_kmers(kmers, file_list_element.filename, tax_id, kmer_len);
		{
			auto seconds_past = std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count();
			if (seconds_past < 1)
				seconds_past = 1;

			size_t megs = total_size/1000000;
			LOG("processed size " << megs << "M = " << (total_size/1000)/seconds_past << "K/sec");
		}
	}

	KmerIO::print_kmers(kmers, kmer_len);

	LOG("total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count());
}


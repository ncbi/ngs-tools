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
#include "kmers.h"
#include "kmer_io.h"
#include "kmer_hash.h"
#include "ready_seq.h"
#include "filename_meta.h"
#include "config_build_index.h"
#include "file_list_loader.h"
#include "build_index.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.37";

size_t weight(size_t kmers_count)
{
	return kmers_count*16;
}

int calculate_window_size_(size_t filesize, bool eukaryota, bool virus)
{
	if (virus)
		return 64;

	if (eukaryota)
		return 8000;

	return 2000;
}

int calculate_window_size(size_t filesize, bool eukaryota, bool virus, int window_divider, int min_window_size)
{
	return std::max(min_window_size, calculate_window_size_(filesize, eukaryota, virus) / window_divider);
}

int main(int argc, char const *argv[])
{
	LOG("build_index version " << VERSION);
	ConfigBuildIndex config(argc, argv);
	LOG("window divider: " << config.window_divider);
	LOG("kmer len: " << config.kmer_len);
	LOG("min window size: " << config.min_window_size);
	LOG("min kmers per clean string: " << config.min_kmers_per_seq);

	auto before = high_resolution_clock::now();

	FileListLoader file_list(config.file_list);

	TaxIdTree tax_id_tree;
	TaxIdTreeLoader::load_tax_id_tree(tax_id_tree, config.tax_parents_file);

	Kmers kmers(tax_id_tree);
	size_t total_size = 0;
	for (auto &file_list_element : file_list.files)
	{
		auto window_size = calculate_window_size(file_list_element.filesize, FilenameMeta::is_eukaryota(file_list_element.filename), FilenameMeta::is_virus(file_list_element.filename), config.window_divider, config.min_window_size);
		auto tax_id = FilenameMeta::tax_id_from(file_list_element.filename);
		LOG(file_list_element.filesize << "\t" << window_size << "\t" << tax_id << "\t" << file_list_element.filename);
		std::vector<std::string> summary;
		total_size += BuildIndex::add_kmers(file_list_element.filename, window_size, config.kmer_len, summary, [&](hash_t kmer){ kmers.add_kmer(kmer, tax_id); });
		{
			auto seconds_past = std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count();
			if (seconds_past < 1)
				seconds_past = 1;

			size_t megs = total_size/1000000;
			LOG("processed size " << megs << "M = " << (total_size/1000)/seconds_past << "K/sec, kmers: " << kmers.storage.size()/1000 << "K, compression rate " << total_size/std::max(size_t(1), weight(kmers.storage.size())));
		}
	}

	KmerIO::print_kmers(kmers, config.kmer_len);

	LOG("total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count());
}

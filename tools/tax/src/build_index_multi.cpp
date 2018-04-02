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
#include "kmers_multi.h"
#include "kmer_io.h"
#include "kmer_hash.h"
#include "ready_seq.h"
#include "config_build_index_multi.h"
#include "file_list_loader.h"
#include "filename_meta.h"
#include "build_index.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.10";

int main(int argc, char const *argv[])
{
	LOG("build_index_multi version " << VERSION);
	Config config(argc, argv);

	auto before = high_resolution_clock::now();

	FileListLoader file_list(config.file_list);

	KmersMulti kmers;
	size_t total_size = 0;
	for (auto &file_list_element : file_list.files)
    {
		auto tax_id = FilenameMeta::tax_id_from(file_list_element.filename);
		LOG(file_list_element.filesize << "\t" << config.window_size << "\t" << tax_id << "\t" << file_list_element.filename);
		total_size += BuildIndex::add_kmers(file_list_element.filename, config.window_size, config.kmer_len, [&](hash_t kmer){ kmers.add_kmer(kmer, tax_id); });
    }

	KmerIO::print_kmers(kmers, config.kmer_len);

	LOG("total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count());
}

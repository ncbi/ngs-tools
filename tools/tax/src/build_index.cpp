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

#include "build_index.h"

const string VERSION = "0.35";

size_t weight(size_t kmers_count)
{
	return kmers_count*16;
}

bool is_eukariota(const string &filename)
{
	return filename.find("/Eukaryota/") != string::npos;
}

bool is_virus(const string &filename)
{
	return filename.find("/Viruses/") != string::npos;
}

int main(int argc, char const *argv[])
{
	cerr << "build_index version " << VERSION << endl;
	ConfigBuildIndex config(argc, argv);
	cerr << "window divider: " << config.window_divider << endl;
	cerr << "kmer len: " << config.kmer_len << endl;

	auto before = high_resolution_clock::now();

	FileListLoader file_list(config.file_list);

	TaxIdTree tax_id_tree;
	TaxIdTreeLoader::load_tax_id_tree(tax_id_tree, config.tax_parents_file);

	Kmers kmers(tax_id_tree);
	size_t total_size = 0;
	for (auto &file_list_element : file_list.files)
	{
		auto window_size = calculate_window_size(file_list_element.filesize, is_eukariota(file_list_element.filename), is_virus(file_list_element.filename), config.window_divider);
		auto tax_id = tax_id_from(file_list_element.filename);
		cerr << file_list_element.filesize << "\t" << window_size << "\t" << tax_id << "\t" << file_list_element.filename << endl;
		total_size += add_kmers(kmers, file_list_element.filename, tax_id, window_size, config.kmer_len);
		{
			auto seconds_past = std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count();
			if (seconds_past < 1)
				seconds_past = 1;

			size_t megs = total_size/1000000;
			cerr << "processed size " << megs << "M = " << (total_size/1000)/seconds_past << "K/sec, kmers: " << kmers.storage.size()/1000 << "K, compression rate " << total_size/std::max(size_t(1), weight(kmers.storage.size())) << endl;
		}
	}

	print_kmers(kmers, config.kmer_len);

	cerr << "total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count() << endl;
}

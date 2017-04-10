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

#include "check_index.h"
#include "tax_id_tree.h"
#include "config_check_index.h"

int main(int argc, char const *argv[])
{
    try
    {
		ConfigCheckIndex config(argc, argv);
		cerr << "check_index version 0.10 " << endl;

		auto before = high_resolution_clock::now();

		FileListLoader file_list(config.file_list);

		TaxIdTree tax_id_tree;
		TaxIdTreeLoader::load_tax_id_tree(tax_id_tree, config.tax_parents_file);

		Kmers kmers(tax_id_tree);
		int kmer_len = load_kmers(kmers, config.kmers_file);
		cerr << "kmer len: " << kmer_len << endl;
		cerr << kmers.storage.size() << " kmers loaded" << endl;

		size_t total_size = 0;
		for (auto &file_list_element : file_list.files)
		{
			auto tax_id = tax_id_from(file_list_element.filename);
			cerr << file_list_element.filesize << "\t" << tax_id << "\t" << file_list_element.filename << endl;
			total_size += check_kmers(kmers, file_list_element.filename, tax_id, kmer_len);
			{
				auto seconds_past = std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count();
				if (seconds_past < 1)
					seconds_past = 1;

				size_t megs = total_size/1000000;
				cerr << "processed size " << megs << "M = " << (total_size/1000)/seconds_past << "K/sec" << endl;
			}
		}

		print_kmers(kmers, kmer_len);

		cerr << "total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count() << endl;

		exit(0); // dont want to wait for KmerMaps destructors
        return 0;
    }
    catch ( const exception & x )
    {
        cerr << x.what() << endl;
//		cerr << "exit 3" << endl;
		return 3;
    }
    catch ( string x )
    {
        cerr << x << endl;
//		cerr << "exit 4" << endl;
		return 4;
    }
    catch ( const char * x )
    {
        cerr << x << endl;
//		cerr << "exit 5" << endl;
		return 5;
    }
    catch ( ... )
    {
        cerr << "unknown exception" << endl;
//		cerr << "exit 6" << endl;
		return 6;
    }
}

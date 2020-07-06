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

#ifndef CONFIG_BUILD_INDEX_H_INCLUDED
#define CONFIG_BUILD_INDEX_H_INCLUDED

#include <string>
#include "log.h"

struct ConfigBuildIndex
{
	std::string file_list, tax_parents_file;
	unsigned int window_divider, kmer_len, min_window_size, min_kmers_per_seq;

	ConfigBuildIndex(int argc, char const *argv[])
	{
		if (argc != 7)
		{
			print_usage();
			exit(1);
		}

		file_list = std::string(argv[1]);
		tax_parents_file = std::string(argv[2]);
		window_divider = std::stoi(std::string(argv[3]));
		kmer_len = std::stoi(std::string(argv[4]));
        min_window_size = std::stoi(std::string(argv[5]));
        min_kmers_per_seq = std::stoi(std::string(argv[6]));
	}

	static void print_usage()
	{
        LOG("need <files.list> <tax.parents> <window divider> <kmer len> <min window size> <min kmers per seq>");
	}
};

#endif

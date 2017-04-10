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

#ifndef CONFIG_CHECK_INDEX_H_INCLUDED
#define CONFIG_CHECK_INDEX_H_INCLUDED

#include <string>

struct ConfigCheckIndex
{
	std::string file_list, tax_parents_file, kmers_file;

	ConfigCheckIndex(int argc, char const *argv[])
	{
		if (argc != 4)
		{
			print_usage();
			exit(1);
		}

		file_list = std::string(argv[1]);
		tax_parents_file = std::string(argv[2]);
		kmers_file = std::string(argv[3]);
	}

	static void print_usage()
	{
		std::cerr << "need <files.list> <tax.parents> <kmers file>" << std::endl;
	}
};

#endif

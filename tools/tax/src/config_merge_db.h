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
#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string file_list, in_ext, out_file;
    size_t expected_size;
	int argc;
	char const **argv;

	std::string arg(int index) const
	{
		if (index >= argc)
			fail();

		return std::string(argv[index]);
	}

	Config(int argc, char const *argv[]) : argc(argc), argv(argv)
	{
		file_list = arg(1);
		in_ext = arg(2);
		if (in_ext == "-")
			in_ext = "";

		out_file = arg(3);
		expected_size = std::stoll(arg(4));
	}

	void fail() const
	{
		print_usage();
        exit(1);
	}

	static void print_usage()
	{
		std::cerr << "need <files.list> <in ext or -> <out file> <expected kmers>" << std::endl;
	}

};


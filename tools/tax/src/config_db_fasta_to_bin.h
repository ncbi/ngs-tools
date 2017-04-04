#ifndef CONFIG_ALIGN_TO_BUILD_INDEX_DB_H_INCLUDED
#define CONFIG_ALIGN_TO_BUILD_INDEX_DB_H_INCLUDED

#include <string>
#include <iostream>

struct Config
{
	std::string fasta_db, out_file;

	Config(int argc, char const *argv[])
	{
		if (argc < 3)
		{
			print_usage();
			throw "";
		}

		fasta_db = argv[1];
		out_file = argv[2];
	}

	static void print_usage()
	{
		std::cerr << "need <fasta db> <out file>" << std::endl;
	}

};

#endif

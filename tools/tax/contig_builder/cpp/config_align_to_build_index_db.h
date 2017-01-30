#ifndef CONFIG_ALIGN_TO_BUILD_INDEX_DB_H_INCLUDED
#define CONFIG_ALIGN_TO_BUILD_INDEX_DB_H_INCLUDED

#include <string>
#include <iostream>
#include <stdexcept>

struct Config
{
	std::string reference, out_file;

	Config(int argc, char const *argv[])
	{
		if (argc < 3)
		{
			print_usage();
			throw std::runtime_error("too few arguments");
		}

		reference = argv[1];
		out_file = argv[2];
	}

	static void print_usage()
	{
		std::cerr << "need <reference fasta> <out file>" << std::endl;
	}

};

#endif

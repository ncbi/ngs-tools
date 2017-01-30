#ifndef CONFIG_KMERS_STAT_H_INCLUDED
#define CONFIG_KMERS_STAT_H_INCLUDED

#include <string>
#include <iostream>
#include <stdexcept>

struct Config
{
	std::string reference;

	Config(int argc, char const *argv[])
	{
		if (argc < 2)
		{
			print_usage();
			throw std::runtime_error("too few arguments");
		}

		reference = argv[1];
	}

	static void print_usage()
	{
		std::cerr << "need <reference fasta>" << std::endl;
	}

};

#endif

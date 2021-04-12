#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string file_list, postfix, allowed_kmers_postfix;
	int argc;
	char const **argv;

	std::string arg(int index) const
	{
		if (index >= argc)
			fail();

		return std::string(argv[index]);
	}

	Config(int argc, char const *argv[]) : argc(argc), argv(argv) // todo: make it right
	{
		file_list = arg(1);
		postfix = arg(2);
		allowed_kmers_postfix = arg(3);
	}

	void fail() const
	{
		print_usage();
        exit(1);
	}

	static void print_usage()
	{
		std::cerr << "need <files.list> <postfix> <allowed kmer postfix>" << std::endl;
	}

};

#endif

#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string file_list;
    int kmer_len;
	int argc;
	char const **argv;

	std::string arg(int index) const
	{
		if (index >= argc)
			fail();

		return std::string(argv[index]);
	}

	Config(int argc, char const *argv[]) : argc(argc), argv(argv), kmer_len(0) // todo: make it right
	{
//		filename = arg(1);
		file_list = arg(1);
		kmer_len = std::stoi(std::string(argv[2]));
	}

	void fail() const
	{
		print_usage();
		throw "";
	}

	static void print_usage()
	{
		std::cerr << "need <files.list> <kmer_len>" << std::endl;
	}

};

#endif

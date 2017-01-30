#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string file_list, tree_file;
    int kmer_len, data_frame, levels_to_save;
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
        tree_file = arg(2);
		kmer_len = std::stoi(arg(3));
		data_frame = std::stoi(arg(4));
        levels_to_save = std::stoi(arg(5));
	}

	void fail() const
	{
		print_usage();
		throw "";
	}

	static void print_usage()
	{
		std::cerr << "need <files.list> <tree file> <kmer_len> <data frame> <levels to save>" << std::endl;
	}

};

#endif

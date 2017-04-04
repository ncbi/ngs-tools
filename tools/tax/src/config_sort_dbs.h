#ifndef CONFIG_SORT_DBS_H_INCLUDED
#define CONFIG_SORT_DBS_H_INCLUDED

#include <string>
#include <iostream>

struct Config
{
	std::string input_filename, out_filename;

	Config(int argc, char const *argv[])
	{
		if (argc < 3)
		{
			print_usage();
			throw "";
		}

		input_filename = argv[1];
		out_filename = argv[2];
	}

	static void print_usage()
	{
		std::cerr << "need <dbs file> <out file>" << std::endl;
	}

};

#endif

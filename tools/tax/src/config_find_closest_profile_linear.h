#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string file_list, profile_file;
	int top_count;
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
		profile_file = arg(2);
		top_count = std::stoi(arg(3));
	}

	void fail() const
	{
		print_usage();
		exit(1);
	}

	static void print_usage()
	{
		std::cerr << "need <files.list> <profile file> <top count>" << std::endl;
	}

};

#endif

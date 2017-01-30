#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string files_list;
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
		files_list = arg(1);
	}

	void fail() const
	{
		print_usage();
		throw "";
	}

	static void print_usage()
	{
		std::cerr << "need <files.list>" << std::endl;
	}

};

#endif

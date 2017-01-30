#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string tree_3mer_file, file_a, file_b;
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
		tree_3mer_file = arg(1);
        file_a = arg(2);
        file_b = arg(3);
	}

	void fail() const
	{
		print_usage();
		throw "";
	}

	static void print_usage()
	{
		std::cerr << "need <tree 3mer file> <file a> <file b>" << std::endl;
	}

};

#endif

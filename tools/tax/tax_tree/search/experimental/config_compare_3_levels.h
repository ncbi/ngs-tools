#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string tree_3mer_file, tree_9mer_file, tree_27mer_file, file_a, file_b;
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
		tree_9mer_file = arg(2);
		tree_27mer_file = arg(3);
        file_a = arg(4);
        file_b = arg(5);
	}

	void fail() const
	{
		print_usage();
		throw "";
	}

	static void print_usage()
	{
		std::cerr << "need <tree 3mer file> <tree 9mer file> <tree 27mer file> <file a> <file b>" << std::endl;
	}

};

#endif

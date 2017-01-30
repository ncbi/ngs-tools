#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string file_list, tree_file, tree_3mer_file;
    int data_frame, levels_to_save;
    int window_len, between_windows_step;
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
        tree_file = arg(2);
        tree_3mer_file = arg(3);

		data_frame = std::stoi(arg(4));
		levels_to_save = std::stoi(arg(5));
        window_len = std::stoi(arg(6));
        between_windows_step = std::stoi(arg(7));
	}

	void fail() const
	{
		print_usage();
		throw "";
	}

	static void print_usage()
	{
		std::cerr << "need <files.list> <tree load/save file> <tree 3mer file> <data frame> <levels to save> <window len> <between_windows_step>" << std::endl;
	}

};

#endif

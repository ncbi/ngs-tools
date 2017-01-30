#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string filename, freq_filename, freq_file_list;
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
		filename = arg(1);
		freq_filename = arg(2);
        if (ends_with(freq_filename, ".list") || freq_filename.find(".list.") != std::string::npos)
        {
            freq_file_list = freq_filename;
            freq_filename = "";
        }
	}

    static bool ends_with(const std::string &s, const std::string &with)
    {
        if (with.size() > s.size()) 
            return false;

        return std::equal(with.rbegin(), with.rend(), s.rbegin());
    }

	void fail() const
	{
		print_usage();
		throw "";
	}

	static void print_usage()
	{
		std::cerr << "need <fasta> <freq_filename>" << std::endl;
	}

};

#endif

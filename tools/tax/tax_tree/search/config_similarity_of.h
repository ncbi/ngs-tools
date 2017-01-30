#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string a_freq_filename, b_freq_filename, freq_file_list;
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
		a_freq_filename = arg(1);
		b_freq_filename = arg(2);
        if (ends_with(b_freq_filename, ".list") || b_freq_filename.find(".list.") != std::string::npos)
        {
            freq_file_list = b_freq_filename;
            b_freq_filename = "";
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
		std::cerr << "need <a freq filename> <b freq filename or files list>" << std::endl;
	}

};

#endif

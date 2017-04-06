#ifndef CONFIG_BUILD_INDEX_H_INCLUDED
#define CONFIG_BUILD_INDEX_H_INCLUDED

#include <string>

struct ConfigBuildIndex
{
	std::string file_list, tax_parents_file;
	unsigned int window_divider, kmer_len;

	ConfigBuildIndex(int argc, char const *argv[])
	{
		if (argc != 5)
		{
			print_usage();
			throw "";
		}

		file_list = std::string(argv[1]);
		tax_parents_file = std::string(argv[2]);
		window_divider = std::stoi(std::string(argv[3]));
		kmer_len = std::stoi(std::string(argv[4]));
	}

	static void print_usage()
	{
		std::cerr << "need <files.list> <tax.parents> <window divider> <kmer len>" << std::endl;
	}
};

#endif

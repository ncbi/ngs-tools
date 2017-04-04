#ifndef CONFIG_CHECK_INDEX_H_INCLUDED
#define CONFIG_CHECK_INDEX_H_INCLUDED

#include <string>

struct ConfigCheckIndex
{
	std::string file_list, tax_parents_file, kmers_file;

	ConfigCheckIndex(int argc, char const *argv[])
	{
		if (argc != 4)
		{
			print_usage();
			throw "";
		}

		file_list = std::string(argv[1]);
		tax_parents_file = std::string(argv[2]);
		kmers_file = std::string(argv[3]);
	}

	static void print_usage()
	{
		std::cerr << "need <files.list> <tax.parents> <kmers file>" << std::endl;
	}
};

#endif

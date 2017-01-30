#ifndef CONFIG_FILTER_DB_H_INCLUDED
#define CONFIG_FILTER_DB_H_INCLUDED

#include <string>

struct ConfigFilterDB
{
	std::string input_file;
	unsigned int only_tax;

	ConfigFilterDB(int argc, char const *argv[]) : only_tax(0)
	{
		if (argc < 2)
		{
			print_usage();
			throw "";
		}

		input_file = std::string(argv[1]);
		if (argc == 4 && std::string(argv[2]) == "-only_tax")
			only_tax = std::stoi(std::string(argv[3]));
	}

	static void print_usage()
	{
		std::cerr << "need <kmers file> [-only_tax <tax_id>]" << std::endl;
	}
};

#endif

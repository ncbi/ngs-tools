#ifndef CONFIG_MIN_COVERAGE_H_INCLUDED
#define CONFIG_MIN_COVERAGE_H_INCLUDED

#include <string>

struct ConfigMinCoverage
{
	std::string reference;
	int kmer_len, lookup_len;

	ConfigMinCoverage(int argc, char const *argv[]) : kmer_len(0), lookup_len(0)
	{
		if (argc != 4)
		{
			print_usage();
			throw "";
		}

		reference = std::string(argv[1]);
		kmer_len = str_to_int(std::string(argv[2]));
		lookup_len = str_to_int(std::string(argv[3]));
	}

	static void print_usage()
	{
		std::cerr << "need <fasta> <kmer_len> <lookup len>" << std::endl;
	}

	static int str_to_int(const std::string &s)
	{
		return std::stoi(s); // throws an exception
	}
};

#endif

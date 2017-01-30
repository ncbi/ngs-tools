#ifndef CONFIG_KMER_DISTRIBUTION_H_INCLUDED
#define CONFIG_KMER_DISTRIBUTION_H_INCLUDED

#include <string>

struct ConfigKmerDistribution
{
	std::string reference;
	int kmer_len, min_coverage;

	ConfigKmerDistribution(int argc, char const *argv[]) : kmer_len(0), min_coverage(0)
	{
		if (argc != 4)
		{
			print_usage();
			throw "";
		}

		reference = std::string(argv[1]);
		kmer_len = str_to_int(std::string(argv[2]));
		min_coverage = str_to_int(std::string(argv[3]));
	}

	static void print_usage()
	{
		std::cerr << "need <text file> <kmer len> <min coverage>" << std::endl;
	}

	static int str_to_int(const std::string &s)
	{
		return std::stoi(s); // throws an exception
	}
};

#endif

#ifndef CONFIG_COMPRESS_H_INCLUDED
#define CONFIG_COMPRESS_H_INCLUDED

#include <string>

struct ConfigCompress
{
	std::string reference, kmer_path_mask;
	int min_dump_coverage;

	ConfigCompress(int argc, char const *argv[]) : min_dump_coverage(0)
	{
		if (argc != 4)
		{
			print_usage();
			throw "";
		}

		reference = std::string(argv[1]);
		kmer_path_mask = std::string(argv[2]);
		min_dump_coverage = str_to_int(std::string(argv[3]));
	}

	static void print_usage()
	{
		std::cerr << "need <text file> <kmers path mask> <min dump coverage>" << std::endl;
	}

	static int str_to_int(const std::string &s)
	{
		return std::stoi(s); // throws an exception
	}
};

#endif

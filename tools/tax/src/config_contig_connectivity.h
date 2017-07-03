#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

#include <string>
#include <iostream>
#include <fstream>
#include <list>

struct Config
{
	std::string fasta_filename, accession;
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
		fasta_filename = arg(1);
		accession = arg(2);
	}

	void fail() const
	{
		print_usage();
		exit(1);
	}

	static void print_usage()
	{
		std::cerr << "need <fasta> <accession>" << std::endl;
	}

};

#endif

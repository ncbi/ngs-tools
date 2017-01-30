#ifndef CONFIG_NEXT_WORDS_H_INCLUDED
#define CONFIG_NEXT_WORDS_H_INCLUDED

#include <string>

struct ConfigNextWords
{
	std::string input;

	ConfigNextWords(int argc, char const *argv[])
	{
		if (argc != 2)
		{
			print_usage();
			throw "";
		}

		input = std::string(argv[1]);
	}

	static void print_usage()
	{
		std::cerr << "need <splitted text file>" << std::endl;
	}
};

#endif

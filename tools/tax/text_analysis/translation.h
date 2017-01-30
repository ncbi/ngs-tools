#ifndef TRANSLATION_H_INCLUDED
#define TRANSLATION_H_INCLUDED

#include <string>
#include <map>

struct Translation
{
	std::map<std::string, char> aminoacid;

	Translation()
	{
		std::string base1 =  "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
		std::string base2 =  "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
		std::string base3 =  "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
		std::string result = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

		for (size_t i =0; i<result.size(); i++)
			aminoacid[std::string(1, base1[i]) + std::string(1, base2[i]) + std::string(1, base3[i])] = result[i];
	}

	void translate(const std::string &nucl_seq, int from, std::string &acids)
	{
		acids.clear();
		acids.reserve((nucl_seq.length() + 2) / 3);
		for (; from + 2 < nucl_seq.length(); from+=3)
			acids.push_back(aminoacid[std::string(&nucl_seq[from], &nucl_seq[from + 3])]);
	}
};

#endif

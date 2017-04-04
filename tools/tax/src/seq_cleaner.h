#ifndef SEQ_CLEANER_H_INCLUDED
#define SEQ_CLEANER_H_INCLUDED

#include <string>
#include <list>
#include "p_string.h"

struct SeqCleaner
{
	typedef std::list<p_string> p_strings;
	p_strings clean_strings;

	bool is_good(char ch)
	{
		return ch == 'A' || ch == 'C' || ch == 'T' || ch == 'G';
	}

	bool is_bad(char ch)
	{
		return !is_good(ch);
	}

	size_t find_good_pos(const std::string &s, size_t from)
	{
		size_t pos = from;
		while ( pos < s.length() && is_bad(s[pos]) )
			pos++;

		return pos;
	}

	size_t find_bad_pos(const std::string &s, size_t from)
	{
		size_t pos = from;
		while ( pos < s.length() && is_good(s[pos]) )
			pos++;

		return pos;
	}

	SeqCleaner(const std::string &s)
	{
		size_t pos = 0;
		while (pos < s.length())
		{
			pos = find_good_pos(s, pos);
			auto bad_pos = find_bad_pos(s, pos);
			auto len = bad_pos - pos;
//			std::cout << pos << " " << bad_pos << std::endl;
			if (len > 0)
			{
				auto p_str = p_string(&s[pos], len);
//				std::cout << std::string(p_str.s) << std::endl;
				clean_strings.push_back(p_str);
			}
			pos = bad_pos;
		}
	}
};


#endif

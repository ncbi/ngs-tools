#ifndef P_STRING_H_INCLUDED
#define P_STRING_H_INCLUDED

#include "stringn.h"

struct p_string
{
	const char *s;
	int len; // todo: unsigned or size_t ?

	p_string(const char *s, int len) : s(s), len(len) 
	{
//			std::cout << s << " " << len << std::endl;
	}

	bool operator < (const p_string &x) const
	{
//			std::cout << x.len << " " << len << std::endl;
		if (x.len != len)
			throw "p_string operator < x.len != len";

		return string_compare(s, x.s, len) < 0;
	}
};


static void print_string(p_string p)
{
	print_string(p.s, p.len);
}


#endif
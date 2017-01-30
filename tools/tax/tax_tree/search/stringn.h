#ifndef STRINGN_H_INCLUDED
#define STRINGN_H_INCLUDED

#include <iostream>

static int string_compare(const char *a, const char *b, unsigned int len)
{
	for (unsigned int i=0; i<len; i++)
	{
		int diff = int(a[i]) - int(b[i]);
		if (diff != 0)
			return diff;
	}

	return 0;
}

static void print_string(const char *s, int len)
{
	for (int i=0; i<len; i++)
		std::cout << s[i];
}

#endif
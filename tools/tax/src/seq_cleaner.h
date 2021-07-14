/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

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

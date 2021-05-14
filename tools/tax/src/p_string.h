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

#ifndef P_STRING_H_INCLUDED
#define P_STRING_H_INCLUDED

#include <ostream>
#include <string>

// todo: update to string_view
struct p_string
{
	const char *s;
	int len;

	p_string(const char *s, int len) : s(s), len(len) {}
    p_string(const std::string &x) : s(x.c_str()), len(x.length()){}

	bool operator < (const p_string &x) const
	{
		if (len != x.len)
			return len < x.len; //throw std::runtime_error("p_string operator < x.len != len");

		return string_compare(s, x.s, len) < 0;
	}

	static int string_compare(const char *a, const char *b, int len)
	{
		for (int i=0; i < len; i++)
		{
			int diff = int(a[i]) - int(b[i]);
			if (diff != 0)
				return diff;
		}

		return 0;
	}
};


std::ostream& operator << (std::ostream &os, const p_string &ps)
{
    for (int i = 0; i < ps.len; i++)
        os << ps.s[i];
    return os;
}
#endif
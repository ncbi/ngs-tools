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
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

#ifndef ACC_LIST_LOADER_H_INCLUDED
#define ACC_LIST_LOADER_H_INCLUDED

#include <string>
#include <fstream>
#include <list>

struct AccListLoader
{
	typedef std::vector<std::string> Files;
	Files files;

	AccListLoader(const std::string &file_list)
	{
		std::ifstream f(file_list);
		if (f.fail() || f.bad())
			throw std::runtime_error(std::string("cannot open file list ") + file_list);

		while (!f.eof())
		{
            std::string s;
            f >> s;
//			std::getline(f, s);
			if (!s.empty())
			    files.push_back(s);
		}
	}
};

#endif

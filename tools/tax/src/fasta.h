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

#ifndef FASTA_H_INCLUDED
#define FASTA_H_INCLUDED

#include <string>
#include <fstream>
#include <vector>

struct Fasta
{
	std::ifstream f;

	std::string last_desc, prev_desc;
	Fasta(const std::string &filename) : f(filename)
	{
		std::string line;
		std::getline(f, line);
		if (line.empty())
			throw "fasta file is empty";

		if (!is_description(line))
			throw "this is not a fasta file";

		last_desc = line;
	}

	static bool is_description(const std::string &s)
	{
		return !s.empty() && s[0] == '>';
	}

	bool get_next_sequence(std::string &result)
	{
		result.clear();
		prev_desc = last_desc;

		std::string line;
		while (!f.eof())
		{
			std::getline(f, line);

			if (line.empty())
				continue;

			if (is_description(line))
			{
				last_desc = line;
				break;
			}

			result += line;
		}

		if (result.empty())
		{
			prev_desc = "";
			return false;
		}

		return true;
	}

	std::string sequence_description() const
	{
		return prev_desc.empty() ? prev_desc : prev_desc.substr(1);
	}

	static size_t filesize(const std::string &filename) 
	{
		std::ifstream f(filename, std::ios_base::binary);
		f.seekg(0, std::ios_base::end);
		return f.tellg();
	}
};

#endif

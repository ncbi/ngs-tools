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

#ifndef TEXT_LOADER_MT_H_INCLUDED
#define TEXT_LOADER_MT_H_INCLUDED

#include <string>
#include <fstream>
#include <stdexcept>

struct TextLoaderSTNoStore
{
	std::ifstream f;

	TextLoaderSTNoStore(const std::string &filename) : f(filename)
	{
		if (f.fail())
			throw std::runtime_error("cannot open input file");
	}

	bool load_next_sequence(std::string &s)
	{
		s.clear();
		while (s.empty() && !f.eof())
			std::getline(f, s);

		return !s.empty();
	}
};

struct FastaWithTaxonomyLoader
{
	std::ifstream f;

	FastaWithTaxonomyLoader(const std::string &filename) : f(filename)
	{
		if (f.fail())
			throw std::runtime_error("cannot open input file");
	}

	bool load_next_sequence(std::string &s, int &tax_id)
	{
		s.clear();
		tax_id = 0;
		if (f.fail() || f.eof())
			return false;

		f >> s;
		f >> tax_id;

		return !s.empty();
	}

};

#endif

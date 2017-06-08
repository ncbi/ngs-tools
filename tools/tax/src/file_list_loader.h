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

#ifndef FILE_LIST_LOADER_H_INCLUDED
#define FILE_LIST_LOADER_H_INCLUDED

#include <string>
#include <fstream>
#include <list>

struct FileListLoader
{
	struct File
	{
		size_t filesize;
		std::string filename;
        File() = default;
        File(size_t filesize, const std::string &filename) : filesize(filesize), filename(filename){}
	};

	typedef std::vector<File> Files;
	Files files;

	FileListLoader(const std::string &file_list)
	{
        if (ends_with(file_list, ".fasta"))
        {
            files.push_back(File(0, file_list));
            return;
        }

		std::ifstream f(file_list);
		if (f.fail() || f.bad())
			throw std::runtime_error(std::string("cannot open file list ") + file_list);

		while (!f.eof())
		{
			File f_rec;
			f >> f_rec.filesize;
			char ch;
			f.get(ch);
			if (ch != '\t')
				throw std::runtime_error("invalid separator in file list - tab expected");

			std::getline(f, f_rec.filename);

			if (f_rec.filename.empty())
				break;

			files.push_back(f_rec);
		}
	}

    static bool ends_with(const std::string &s, const std::string &end)
    {
        if (end.size() > s.size()) 
            return false;

        return std::equal(end.rbegin(), end.rend(), s.rbegin());
    }
};

#endif

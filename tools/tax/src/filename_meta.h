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
#pragma once
#include <string>

struct FilenameMeta
{
    static bool is_eukaryota(const std::string &filename)
    {
	    return filename.find("/Eukaryota/") != std::string::npos;
    }

    static bool is_virus(const std::string &filename)
    {
	    return filename.find("/Viruses/") != std::string::npos;
    }

    static int tax_id_from(const std::string &filename)
    {
	    auto to = filename.find_last_of('.');
	    auto from = filename.find_last_of('/') + 1;
	    return stoi(filename.substr(from, to - from));
    }
/*
    static int tax_id_from_folder(const std::string &filename)
    {
	    auto to = filename.find_last_of('/');
	    auto from = filename.find_last_of('/', to - 1) + 1;
//		std::cout << filename << " " << from << " " << to << " |" << filename.substr(from, to - from) << "|" << std::endl;
	    return stoi(filename.substr(from, to - from));
    }
*/

};

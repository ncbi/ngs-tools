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

#ifndef DBS_H_INCLUDED
#define DBS_H_INCLUDED

#include "database_io.h"
#include <string>
#include <fstream>
#include <iostream>

struct DBS
{
	#pragma pack(push)
	#pragma pack(4)

	struct KmerTax
	{
		hash_t kmer;
		int tax_id;

		KmerTax(hash_t kmer = 0, int tax_id = 0) : kmer(kmer), tax_id(tax_id){}
	};

	#pragma pack(pop)

	typedef std::vector<KmerTax> Kmers;
};

struct DBSIO
{
	static const int VERSION = 1;

	struct DBSHeader
	{
		size_t version, kmer_len;
		DBSHeader(size_t kmer_len = 0) : version(VERSION), kmer_len(kmer_len){}
	};

	template <class C>
	static void save_dbs(const std::string &out_file, const std::vector<C> &kmers, size_t kmer_len)
	{
		std::ofstream f(out_file);
		DBSHeader header(kmer_len);
        IO::write(f, header);
        IO::save_vector(f, kmers);
	}

	template <class C>
	static size_t load_dbs(const std::string &filename, std::vector<C> &kmers)
	{
    	std::ifstream f(filename, std::ios::binary | std::ios::in);
	    if (f.fail() || f.eof())
		    throw std::runtime_error(std::string("cannot load dbs ") + filename);

        DBSHeader header;
        IO::read(f, header);
		if (header.version != VERSION)
			throw std::runtime_error("unsupported dbs file version");

        IO::load_vector(f, kmers);

		if (header.kmer_len < 1 || header.kmer_len > 64)
			throw std::runtime_error("load_dbs:: invalid kmer_len");

		return header.kmer_len;
	}
};

#endif

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
	static const int VERSION = 1;
//	typedef uint64_t hash_t;

	#pragma pack(push)
	#pragma pack(4)

	struct KmerTax
	{
		hash_t kmer;
		int tax_id;

		KmerTax(hash_t kmer = 0, int tax_id = 0) : kmer(kmer), tax_id(tax_id){}

		//bool operator < (const KmerTax &x) const
		//{
		//	if (tax_id == x.tax_id)
		//		return kmer < x.kmer;

		//	return tax_id < x.tax_id;
		//}
	};

	#pragma pack(pop)

	typedef std::vector<KmerTax> Kmers;
//	typedef vector<hash_t> Hashes;

	struct DBSHeader
	{
		size_t version, kmer_len;
		DBSHeader(size_t kmer_len = 0) : version(VERSION), kmer_len(kmer_len){}
	};

	template <class C>
	static void save_dbs(const std::string &out_file, const std::vector<C> &kmers, size_t kmer_len)
	{
		std::ofstream f(out_file);
		f.flush();

		size_t kmers_size = kmers.size();
//		cerr << "kmer len: " << kmer_len << endl;
//		cerr << "kmers: " << kmers_size << endl;
		DBSHeader header(kmer_len);

		f.write((char*)&header, sizeof(header));
		f.write((char*)&kmers_size, sizeof(kmers_size));
		f.write((char*)&kmers[0], kmers.size()*sizeof(kmers[0]));
	}

	template <class C>
	static size_t load_dbs(const std::string &filename, std::vector<C> &kmers)
	{
		DBSHeader header;
		load_structure(filename, header);
		std::cerr << "version: " << header.version << std::endl;
		std::cerr << "kmer len: " << header.kmer_len << std::endl;

		if (header.version != VERSION)
			throw std::runtime_error("unsupported version");

		if (header.kmer_len < 10 || header.kmer_len > 64)
			throw std::runtime_error("invalid kmer_len");

		load_vector(filename, kmers, sizeof(header));
		return header.kmer_len;
	}
};

#endif

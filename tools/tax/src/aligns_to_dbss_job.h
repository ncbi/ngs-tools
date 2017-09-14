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

#ifndef ALIGNS_TO_DBSS_JOB_H_INCLUDED
#define ALIGNS_TO_DBSS_JOB_H_INCLUDED

#include "aligns_to_dbs_job.h"
#include <set>

struct DBSSJob : public DBSJob
{
	DBSSJob(const Config &config) : DBSJob(config)
	{
		DBSIO::DBSHeader header;

	    std::ifstream f(config.dbss, std::ios::binary | std::ios::in);
	    if (f.fail() || f.eof())
		    throw std::runtime_error(std::string("cannot open dbss ") + config.dbss);

		IO::read(f, header);
		kmer_len = header.kmer_len;

		DBSAnnotation annotation;
		auto sum_offset = load_dbs_annotation(config.dbss + ".annotation", annotation);
		if (sum_offset != IO::filesize(config.dbss))
			throw std::runtime_error("inconsistent dbss annotation file");

		auto tax_list = load_tax_list(config.dbss_tax_list);
		if (tax_list.empty())
			throw std::runtime_error("empty tax list");

		load_dbss(config.dbss, tax_list, annotation);
	}

	typedef int tax_id_t;

	struct DBSAnnot
	{
		tax_id_t tax_id;
		size_t count, offset;

		DBSAnnot(int tax_id, size_t count, size_t offset) : tax_id(tax_id), count(count), offset(offset){}

		bool operator < (const DBSAnnot &x) const
		{
			return tax_id < x.tax_id;
		}
	};

	typedef std::vector<DBSAnnot> DBSAnnotation;
	typedef std::vector<tax_id_t> TaxList;

	static size_t load_dbs_annotation(const std::string &filename, DBSAnnotation &annotation)
	{
		std::ifstream f(filename);
		if (f.fail())
			throw std::runtime_error("cannot open annotation file");

		size_t offset = sizeof(DBSIO::DBSHeader) + sizeof(size_t);
		tax_id_t prev_tax = 0;

		while (!f.eof())
		{
			DBSAnnot a(0, 0, 0);
			f >> a.tax_id >> a.count;
			if (f.fail())
				break;

			if (!a.count)
				throw std::runtime_error("bad annotation format - bad count");

			if (prev_tax >= a.tax_id)
				throw std::runtime_error("bad annotation format - tax less");

			a.offset = offset;
			annotation.push_back(a);
			offset += sizeof(hash_t) *a.count;
			prev_tax = a.tax_id;
		}

		return offset;
	}

	static TaxList load_tax_list(const std::string &filename)
	{
		std::ifstream f(filename);
		if (f.fail())
			throw std::runtime_error("cannot open tax list file");

		TaxList taxes;

		while (!f.eof())
		{
			tax_id_t t = 0;
			f >> t;
			if (f.fail())
				break;

			taxes.push_back(t);
		}

		sort(taxes.begin(), taxes.end());
		return taxes;
	}

    void load_dbss(const std::string &filename, const TaxList &tax_list, const DBSAnnotation &annotation)
	{
        assert(hash_array.empty());

        std::ifstream f(filename);
        if (f.fail() || f.eof())
            throw std::runtime_error("cannot open dbss file");
        
        size_t total_hashes_count = 0;
        for (auto tax_id : tax_list) {
            for (auto& annot : annotation) {
                if (annot.tax_id == tax_id) {
                    total_hashes_count += annot.count;
                }
            }
        }
        hash_array.reserve(total_hashes_count);
        
        for (auto tax_id : tax_list) {
            for (auto& annot : annotation) {
                if (annot.tax_id == tax_id && annot.count > 0) {
                    std::vector<hash_t> hashes;
                    IO::load_vector_no_size(f, hashes, annot.offset, annot.count);
                    for (auto hash : hashes) {
                        hash_array.emplace_back(hash, annot.tax_id);
                    }
                }
            }
        }
        assert(hash_array.size() == total_hashes_count);
        
        LOG("dbss parts loaded (" << (total_hashes_count / 1000 / 1000) << "m kmers)");
        assert(!hash_array.empty());
        std::sort(hash_array.begin(), hash_array.end());
        LOG("dbss parts merged");
	}
};

#endif

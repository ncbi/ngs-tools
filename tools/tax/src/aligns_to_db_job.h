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

#ifndef ALIGNS_TO_DB_JOB_H_INCLUDED
#define ALIGNS_TO_DB_JOB_H_INCLUDED

#include "aligns_to_job.h"
#include "hash.h"
#include "seq_transform.h"
#include <map>
#include "omp_adapter.h"

struct DBJob : public Job
{
	typedef std::vector<hash_t> HashSortedArray;

	HashSortedArray hash_array;
	size_t kmer_len;
	const Config &config;

	DBJob(const Config &config) : config(config)
	{
		kmer_len = DBSIO::load_dbs(config.db, hash_array);
	}

	struct Matcher
	{
		const HashSortedArray &hash_array;
		size_t kmer_len;
		Matcher(const HashSortedArray &hash_array, size_t kmer_len) : hash_array(hash_array), kmer_len(kmer_len){}

		int operator() (const std::string &seq) const 
		{
			int found = 0;
			Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
				{
					if (in_db(hash) > 0)
						found++;

					return !found;
				});

			return found;
		}

		bool in_db(hash_t hash) const
		{
			hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
			return std::binary_search(hash_array.begin(), hash_array.end(), hash);
		}
	};

	virtual void run(const std::string &filename, std::ostream &out_f)
	{
		Matcher m(hash_array, kmer_len);
		BasicPrinter print(out_f);
		Job::run<Matcher, BasicPrinter>(filename, print, m, kmer_len, config.spot_filter_file, config.unaligned_only);
	}
};

#endif

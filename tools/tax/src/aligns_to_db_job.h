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
#include <list>

struct BasicMatchId
{
	int seq_id;

	BasicMatchId(int seq_id, int matches) : seq_id((int)seq_id){}
	bool operator < (const BasicMatchId &b) const { return seq_id < b.seq_id; }
};

struct KmerBasicMatchId
{
    struct Matches : public std::list<hash_t>
    {
        operator bool() const { return !empty(); }
    };

	int seq_id;
    const Matches matches;

	KmerBasicMatchId(int seq_id, const Matches &matches) : seq_id((int)seq_id), matches(matches) {}
	bool operator < (const KmerBasicMatchId &b) const { return seq_id < b.seq_id; }
};

struct BasicPrinter
{
    IO::Writer &writer;
    BasicPrinter(IO::Writer &writer) : writer(writer){}

	void operator() (const std::vector<Reader::Fragment> &processing_sequences, const std::vector<BasicMatchId> &ids)
	{
		for (auto seq_id : ids)
			writer.f() << processing_sequences[seq_id.seq_id].spotid << std::endl;

        writer.check();
	}
};

struct KmerBasicPrinter
{
    IO::Writer &writer;
    int kmer_len = 0;
    KmerBasicPrinter(IO::Writer &writer, int kmer_len) : writer(writer), kmer_len(kmer_len){}

	void operator() (const std::vector<Reader::Fragment> &processing_sequences, const std::vector<KmerBasicMatchId> &ids)
	{
		for (auto seq : ids)
		for (auto kmer : seq.matches)
		    writer.f() << Hash<hash_t>::str_from_hash(kmer, kmer_len) << std::endl;

        writer.check();
	}
};

struct DBJob : public Job
{
	typedef std::vector<hash_t> HashSortedArray;

	HashSortedArray hash_array;
	size_t kmer_len;

	DBJob(const std::string &db)
	{
		kmer_len = DBSIO::load_dbs(db, hash_array);
	}

	struct Matcher
	{
		const HashSortedArray &hash_array;
		size_t kmer_len;
		Matcher(const HashSortedArray &hash_array, size_t kmer_len) : hash_array(hash_array), kmer_len(kmer_len){}

		int operator() (const std::string &seq) const 
		{
			int found = 0;
			Hash<hash_t>::for_all_hashes_do(seq, (int)kmer_len, [&](hash_t hash)
				{
					if (in_db(hash))
						found++;

					return !found;
				});

			return found;
		}

		bool in_db(hash_t hash) const
		{
			hash = seq_transform<hash_t>::min_hash_variant(hash, (int)kmer_len);
			return std::binary_search(hash_array.begin(), hash_array.end(), hash);
		}
	};

	struct KmerMatcher
	{
		const HashSortedArray &hash_array;
		size_t kmer_len;
		KmerMatcher(const HashSortedArray &hash_array, size_t kmer_len) : hash_array(hash_array), kmer_len(kmer_len){}

		KmerBasicMatchId::Matches operator() (const std::string &seq) const 
		{
            KmerBasicMatchId::Matches matches; // todo: optimize. though not really urgent
			Hash<hash_t>::for_all_hashes_do(seq, (int)kmer_len, [&](hash_t hash)
				{
					if (in_db(hash))
						matches.push_back(hash);

					return true;
				});

			return matches;
		}

		bool in_db(hash_t hash) const
		{
			hash = seq_transform<hash_t>::min_hash_variant(hash, (int)kmer_len);
			return std::binary_search(hash_array.begin(), hash_array.end(), hash);
		}
	};

    bool print_kmers_only = false;

	virtual void run(const std::string &filename, IO::Writer &writer, const Config &config) override
	{
        print_kmers_only = config.print_kmers_only;
		Job::run_for_matcher(filename, config.spot_filter_file, config.unaligned_only, config.optimization_ultrafast_skip_reader, config.chunk_size, [&](const std::vector<Reader::Fragment> &chunk){ match_and_print_chunk(chunk, writer); } );
	}

    virtual void match_and_print_chunk(const std::vector<Reader::Fragment> &chunk, IO::Writer &writer)
    {
        if (print_kmers_only)
        {
		    KmerMatcher matcher(hash_array, kmer_len);
    		KmerBasicPrinter print(writer, kmer_len);
            Job::match_and_print<KmerMatcher, KmerBasicPrinter, KmerBasicMatchId>(chunk, print, matcher);
        }
        else
        {
    		Matcher matcher(hash_array, kmer_len);
    		BasicPrinter print(writer);
            Job::match_and_print<Matcher, BasicPrinter, BasicMatchId>(chunk, print, matcher);
        }
    }
};

#endif

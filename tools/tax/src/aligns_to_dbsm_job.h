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

#include "aligns_to_dbs_job.h"

struct DBSMJob : public Job
{
	struct KmerTax : public DBS::KmerTaxMulti
	{
        KmerTax() = default;
		KmerTax(hash_t kmer, const std::vector<int> &tax_ids) : DBS::KmerTaxMulti(kmer, tax_ids) { } // todo: remove constructor from KmerTax for faster loading ?

		bool operator < (const KmerTax &x) const // for binary search by hash
		{
			return kmer < x.kmer;
		}
	};

	typedef std::vector<KmerTax> HashSortedArray;

	HashSortedArray hash_array;
	size_t kmer_len = 0;

public:

	DBSMJob(const std::string &dbsm)
	{
		kmer_len = DBSIO::load_dbsm(dbsm, hash_array);
	}

    typedef DBSJob::Hits Hits;

	virtual size_t db_kmers() const override { return hash_array.size();}

	struct Matcher
	{
		const HashSortedArray &hash_array;
		int kmer_len;
        const std::vector<int> EMPTY_TAXES;

		Matcher(const HashSortedArray &hash_array, int kmer_len) : hash_array(hash_array), kmer_len(kmer_len) { }

        const std::vector<int> &find_hash(hash_t hash, const std::vector<int> &default_value) const
        {
            auto first = hash_array.begin();
            auto last = hash_array.end();
            first = std::lower_bound(first, last, KmerTax(hash, default_value));
            return ((first == last) || (hash < first->kmer) ) ? default_value : first->tax_ids;
        }

		Hits operator() (const std::string &seq) const 
		{
			Hits hits;
			Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
				{
					auto &tax_ids = get_db_tax(hash);
                    for (auto tax_id : tax_ids)
                        hits[tax_id] ++;

					return true;
				});

			return hits;
		}

		const std::vector<int> &get_db_tax(hash_t hash) const
		{
			hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
            return find_hash(hash, EMPTY_TAXES);
		}
	};

    typedef DBSJob::TaxMatchId TaxMatchId;
    typedef DBSJob::TaxPrinter TaxPrinter;

    bool hide_counts = false;

	virtual void run(const std::string &filename, IO::Writer &writer, const Config &config) override
	{
        hide_counts = config.hide_counts;
		Job::run_for_matcher(filename, config.spot_filter_file, config.unaligned_only, config.optimization_ultrafast_skip_reader, [&](const std::vector<Reader::Fragment> &chunk) { match_and_print_chunk(chunk, writer); } );
	}

    virtual void match_and_print_chunk(const std::vector<Reader::Fragment> &chunk, IO::Writer &writer)
    {
		Matcher m(hash_array, (int)kmer_len);
		TaxPrinter print(!hide_counts, false, writer);
		Job::match_and_print<Matcher, TaxPrinter, TaxMatchId>(chunk, print, m);
    }
};


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

#ifndef ALIGNS_TO_DBS_JOB_H_INCLUDED
#define ALIGNS_TO_DBS_JOB_H_INCLUDED

struct DBSJob : public Job
{
//	typedef uint64_t hash_t;
	struct KmerTax : public DBS::KmerTax
	{
		KmerTax(hash_t kmer = 0, int tax_id = 0) : DBS::KmerTax(kmer, tax_id){} // todo: remove constructor from KmerTax for faster loading ?

		bool operator < (const KmerTax &x) const // for binary search by hash
		{
			return kmer < x.kmer;
		}
	};

//	struct KmerTax
//	{
//		hash_t kmer;
//		int tax_id;
////		KmerTax(hash_t kmer, int tax_id) : kmer(kmer), tax_id(tax_id){} // do not want it to be called when loading
//	static KmerTax kmer_tfrom(hash_t kmer, int tax_id) { KmerTax t; t.kmer = kmer; t.tax_id = tax_id; return t; }
//
//	};
//
	const Config &config;
	typedef vector<KmerTax> HashSortedArray;

	HashSortedArray hash_array;
    static const int DEFAULT_KMER_LEN = 32;
	typedef unsigned int tax_t;
	size_t kmer_len;

protected:
	DBSJob(const Config &config) : kmer_len(0), config(config)
	{
	}

public:

	struct Hits : public map<tax_t, int>
	{
		operator bool() const { return size() != 0; }

        void operator += (const Hits &x)
        {
            for (auto &other : x)
                (*this)[other.first] += other.second;
        }
	};

	virtual size_t db_kmers() const { return hash_array.size();}

	struct Matcher
	{
        typedef vector<pair<size_t, size_t> > HashLookupTable;
		const HashSortedArray &hash_array;
        HashLookupTable hash_lookup_table;
        int hash_lookup_shift;
		int kmer_len;
		Matcher(const HashSortedArray &hash_array, int kmer_len) : hash_array(hash_array), kmer_len(kmer_len)
        {
            // determining size of lookup key
            int lookup_key_bits = 1;
            while ((hash_array.size() >> lookup_key_bits) > 5)
            {
                lookup_key_bits += 1;
            }
            hash_lookup_shift = kmer_len * 2 - lookup_key_bits;

            const size_t bucket_count = size_t(1) << lookup_key_bits;
            std::cerr << "creating lookup table with " << bucket_count << " buckets, on average " << (float(hash_array.size()) / bucket_count) << " hashes per bucket" << std::endl;
            hash_lookup_table.resize(bucket_count);

            // figuring out bucket ranges
            size_t hash_idx = 0;
            hash_t last_hash = 0;
            for (size_t bucket_idx = 0; bucket_idx < bucket_count; ++bucket_idx)
            {
                auto &bucket = hash_lookup_table[bucket_idx];
                bucket.first = hash_idx;
                while (true)
                {
                    if (hash_idx >= hash_array.size())
                        break;
                    hash_t hash = hash_array[hash_idx].kmer;
                    assert(hash >= last_hash);
                    const size_t hash_bucket = hash >> hash_lookup_shift;
                    if (hash_bucket != bucket_idx)
                        break;
                    
                    ++hash_idx;
                    last_hash = hash;
                }
                bucket.second = hash_idx;
                //std::cerr << "bucket " << bucket_idx << " " << (bucket.second - bucket.first) << " members, range " << bucket.first << " - " << bucket.second << std::endl;
            }
        }

        int find_hash(hash_t hash, int default_value) const
        {
            hash_t bucket_idx = hash >> hash_lookup_shift;
            auto &bucket = hash_lookup_table[bucket_idx];
            auto first = hash_array.begin() + bucket.first;
            auto last = hash_array.begin() + bucket.second;
            first = std::lower_bound(first, last, KmerTax(hash, 0));
            return ((first == last) || (hash < first->kmer) ) ? default_value : first->tax_id;
        }

		Hits operator() (const string &seq) const 
		{
			Hits hits;
			Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
				{
					if (auto tax_id = get_db_tax(hash))
					{
						hits[tax_id] ++;
//						cerr << tax_id << " " << Hash<hash_t>::str_from_hash(hash, KMER_LEN) << endl;
					}

					return true;
				});

			return hits;
		}

		tax_t get_db_tax(hash_t hash) const
		{
			auto tax_id = get_db_tax_0_variations(hash);
#if 0
			if (tax_id)
				return tax_id;
			
			seq_transform<hash_t>::for_all_1_char_variations_do(hash, KMER_LEN, [&](hash_t hash)
				{
					hash = seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN);
					tax_id = find_hash(hash, 0);
					return !tax_id;
				});

			if (tax_id)
				return tax_id;

			seq_transform<hash_t>::for_all_2_char_variations_do(hash, KMER_LEN, [&](hash_t hash)
				{
					hash = seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN);
					tax_id = find_hash(hash, 0);
					return !tax_id;
				});

			if (tax_id)
				return tax_id;

			seq_transform<hash_t>::for_all_3_char_variations_do(hash, KMER_LEN, [&](hash_t hash)
				{
					hash = seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN);
					tax_id = find_hash(hash, 0);
					return !tax_id;
				});
#endif
			return tax_id;
		}

		tax_t get_db_tax_0_variations(hash_t hash) const
		{
			hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
            return find_hash(hash, 0);
		}
	};

	struct TaxMatchId
	{
		int seq_id;
		Hits hits;
		TaxMatchId(int seq_id, const Hits &hits) : seq_id(seq_id), hits(hits)	{}

		bool operator < (const TaxMatchId &b) const { return seq_id < b.seq_id; }
	};

	struct TaxPrinter
	{
		std::ostream &out_f;
        const bool print_counts;
		TaxPrinter(std::ostream &out_f, bool print_counts) : out_f(out_f), print_counts(print_counts) {}

		void operator() (const std::vector<Reader::Fragment> &processing_sequences, const vector<TaxMatchId> &ids)
		{
			for (auto seq_id : ids)
			{
                for (auto c : processing_sequences[seq_id.seq_id].spotid) {
                    if (c == '\t' || c == '\n') {
                        c = ' ';
                    }
                    out_f << c;
                }
                for (auto &hit : seq_id.hits) {
                    out_f << '\t' << hit.first;
                    if (print_counts && hit.second > 1)
                        out_f << 'x' << hit.second;
                }
                out_f << std::endl;
			}
        }
	};

	virtual void run(const string &filename, std::ostream &out_f)
	{
		Matcher m(hash_array, kmer_len);
		TaxPrinter print(out_f, config.print_counts);
		Job::run<Matcher, TaxPrinter, TaxMatchId>(filename, print, m, kmer_len, config.spot_filter_file, config.unaligned_only);
	}
};

struct DBSBasicJob : public DBSJob
{
	DBSBasicJob(const Config &config) : DBSJob(config)
	{
		kmer_len = DBS::load_dbs(config.dbs, hash_array);
	}
};

#endif

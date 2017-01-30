#ifndef ALIGNS_TO_DB_JOB_H_INCLUDED
#define ALIGNS_TO_DB_JOB_H_INCLUDED

struct DBJob : public Job
{
	typedef vector<hash_t> HashSortedArray;

	HashSortedArray hash_array;
	size_t kmer_len;
	const Config &config;

	DBJob(const Config &config) : config(config)
	{
		kmer_len = DBS::load_dbs(config.db, hash_array);
	}

	struct Matcher
	{
		const HashSortedArray &hash_array;
		size_t kmer_len;
		Matcher(const HashSortedArray &hash_array, size_t kmer_len) : hash_array(hash_array), kmer_len(kmer_len){}

		int operator() (const string &seq) const 
		{
			int found = 0;
			Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
				{
					if (in_db(hash) > 0)
						found++;

					return COUNT_MATCHES || !found;
				});

			return found;
		}

		bool in_db(hash_t hash) const
		{
			hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
			return std::binary_search(hash_array.begin(), hash_array.end(), hash);
		}
	};

	virtual void run(const string &filename, std::ostream &out_f)
	{
		Matcher m(hash_array, kmer_len);
		BasicPrinter print(out_f);
		Job::run<Matcher, BasicPrinter>(filename, print, m, kmer_len, config.spot_filter_file, config.unaligned_only);
	}
};

#endif

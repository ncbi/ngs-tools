#ifndef KMER_LOADER_H_INCLUDED
#define KMER_LOADER_H_INCLUDED

#include "mem_usage.h"
#include "reader.h"
#include "vdb_reader.h"

struct KmerLoader
{
	KmerMap32 &kmers;
    Reader::Params reader_params;

	KmerLoader(KmerMap32 &kmers, bool unaligned_only, const std::string& filter_file, bool exclude_filter) : 
		kmers(kmers)
    {
        reader_params.filter_file = filter_file;
        reader_params.exclude_filter = exclude_filter;
        reader_params.unaligned_only = unaligned_only;
        reader_params.read_qualities = false;
    };

	void load(const std::string &accession)
	{
		auto before = std::chrono::high_resolution_clock::now();
        load_32(accession);
		std::cerr << "loading total time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - before ).count() << std::endl;
	}

	template <class hash_t>
	struct NoCheck
	{
		bool operator () (hash_t hash){ return true; }
	};

	void load_32(const std::string &accession)
	{
		auto before = std::chrono::high_resolution_clock::now();
		load_min_mem_map<KmerMap32>(accession, kmers, NoCheck<KmerMap32::hash_t>());
		std::cerr << "32mer loading time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - before ).count() << std::endl;
		std::cerr << "32mer real size: " << kmers.size() << std::endl;

		before = std::chrono::high_resolution_clock::now();
		std::cerr << "mem usage " << mem_usage()/1000000000 << "G" << std::endl;
		kmers.optimize();
		std::cerr << "32mer optimization time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - before ).count() << std::endl;
		std::cerr << "32mer optimized size: " << kmers.size() << std::endl;
		std::cerr << "mem usage " << mem_usage()/1000000000 << "G" << std::endl;
	}

	template <class KmerMap, class Predicate>
	void load_min_mem_map(const std::string &accession, KmerMap &kmers, Predicate pred)
	{
        auto reader = Reader::create(accession, reader_params);

        const int THREADS = 4;
        #pragma omp parallel num_threads(THREADS)
        {
            std::vector<Reader::Fragment> chunk;
            bool done = false;
            while (!done) 
            {
                #pragma omp critical (read)
                {
                    done = !reader->read_many(chunk);
                }

                for (auto& frag: chunk) 
                {
                    auto lambda = [&](typename KmerMap::hash_t hash) 
                    {
                        if (pred(hash))
                            kmers.add(hash);
                        return true;
                    };

                    Hash<typename KmerMap::hash_t>::for_all_hashes_do(frag.bases, kmers.kmer_len, lambda);
                }
            }
        }
	}
};


#endif

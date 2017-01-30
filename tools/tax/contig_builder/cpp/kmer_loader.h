#ifndef KMER_LOADER_H_INCLUDED
#define KMER_LOADER_H_INCLUDED

#include "mem_usage.h"
#include "reader.h"
#include "vdb_reader.h"

struct KmerLoader
{
	KmerMap64 &kmers64;
	KmerMap16 &kmers16;
    Reader::Params reader_params;

	KmerLoader(KmerMap64 &kmers64, KmerMap16 &kmers16, bool unaligned_only, const std::string& filter_file, bool exclude_filter) : 
		kmers64(kmers64), kmers16(kmers16)
    {
        reader_params.filter_file = filter_file;
        reader_params.exclude_filter = exclude_filter;
        reader_params.unaligned_only = unaligned_only;
        reader_params.read_qualities = false;
    };

	void load(const std::string &accession)
	{
		auto before = std::chrono::high_resolution_clock::now();
#if 0
		load_16(accession);
#endif
        load_64(accession);
		std::cerr << "loading total time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - before ).count() << std::endl;
	}

	template <class hash_t>
	struct NoCheck
	{
		bool operator () (hash_t hash){ return true; }
	};

	void load_16(const std::string &accession)
	{
		size_t run_bases = 200 * VdbReader(accession).stats().expected_spot_count; // todo: tune

		auto before = std::chrono::high_resolution_clock::now();
		size_t expected_16_size = (size_t)std::min(size_t(1000000000), run_bases/2);
		std::cerr << "16 expected size: " << expected_16_size << std::endl;

		kmers16.reserve(expected_16_size);
		load_min_mem_map<KmerMap16>(accession, kmers16, NoCheck<KmerMap16::hash_t>());
		std::cerr << "loading 16 time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - before ).count() << std::endl;
		std::cerr << "16 real size: " << kmers16.size() << std::endl;

		before = std::chrono::high_resolution_clock::now();
		kmers16.optimize();
		std::cerr << "16 optimization time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - before ).count() << std::endl;
		std::cerr << "16 optimized size: " << kmers16.size() << std::endl;
		std::cerr << "mem usage " << mem_usage()/1000000000 << "G" << std::endl;
	}

	void load_64(const std::string &accession)
	{
		auto before = std::chrono::high_resolution_clock::now();
//		auto expected_64_size = (kmers32.size()/10)*13;
	//	std::cerr << "64 expected size: " << expected_64_size << std::endl;
//		kmers64.reserve(expected_64_size);
		load_min_mem_map<KmerMap64>(accession, kmers64, NoCheck<KmerMap64::hash_t>());
		std::cerr << "loading 64 time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - before ).count() << std::endl;
		std::cerr << "64 real size: " << kmers64.size() << std::endl;

		before = std::chrono::high_resolution_clock::now();
		std::cerr << "mem usage " << mem_usage()/1000000000 << "G" << std::endl;
		kmers64.optimize();
		std::cerr << "64 optimization time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - before ).count() << std::endl;
		std::cerr << "64 optimized size: " << kmers64.size() << std::endl;
		std::cerr << "mem usage " << mem_usage()/1000000000 << "G" << std::endl;
	}

	template <class KmerMap, class Predicate>
	void load_min_mem_map(const std::string &accession, KmerMap &kmers, Predicate pred)
	{
        auto reader = Reader::create(accession, reader_params);

        #pragma omp parallel
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

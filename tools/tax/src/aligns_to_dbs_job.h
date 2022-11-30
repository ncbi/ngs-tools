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

#include "aligns_to_job.h"
#include <set>
#include <map>
#include "hash.h"
#include "seq_transform.h"
#include "p_string.h"
#include "dbs.h"
#include <mutex>
#include "tax_collator.hpp"

// incompatible with multiple intput files option todo: fix by moving table creation to constructor and enable
#define LOOKUP_TABLE 1

struct DBSJob : public Job
{
    struct KmerTax : public DBS::KmerTax
    {
        KmerTax(hash_t kmer = 0, int tax_id = 0) : DBS::KmerTax(kmer, tax_id) { } // todo: remove constructor from KmerTax for faster loading ?

        bool operator < (const KmerTax &x) const // for binary search by hash
        {
            return kmer < x.kmer;
        }
    };

    typedef std::set<hash_t> hash_set;
    typedef std::vector<KmerTax> HashSortedArray;

    HashSortedArray hash_array;
    typedef unsigned int tax_t; // todo: remove duplicate definition of tax_t and tax_id_t
    size_t kmer_len = 0;

    public:

    struct Hits : public std::map<tax_t, hash_set>
    {
        operator bool() const { return size() != 0; }

        void operator += (const Hits &x)
        {
            for (auto &other : x)
                // (*this)[other.first] += other.second;
                for (auto other_khit : other.second){
                    (*this)[other.first].emplace(other_khit);
                    //(*this)[other.first].insert(other.second.begin(),other.second.end());
                }
        }
        void print_uniq_hits(IO::Writer &uniq_file)
        {
            if (uniq_file.stream_f){
                for (auto tax_hit : *this){
                    uniq_file.f() <<  tax_hit.first << '\t' << tax_hit.second.size() << std::endl;
                }
            }
            else {
                std::cout << "Unique" << unique_hits.size() << std::endl;
                for (auto tax_hit : unique_hits){
                    std::cout <<  tax_hit.first << '\t' << tax_hit.second.size() << std::endl;
                }
            }
        }

    };
    static Hits unique_hits;

    virtual size_t db_kmers() const override { return hash_array.size();}

    struct Matcher
    {
#if LOOKUP_TABLE
        typedef std::vector<size_t> HashLookupTable;
        HashLookupTable hash_lookup_table;
        int hash_lookup_shift;
#endif

        const HashSortedArray &hash_array;
        int kmer_len;
        int max_lookups_per_seq = 0;
        bool unique = false;

        Matcher(const HashSortedArray &hash_array, int kmer_len, int max_lookups_per_seq, bool unique) : hash_array(hash_array), kmer_len(kmer_len), max_lookups_per_seq(max_lookups_per_seq), unique(unique)
        {
            if (max_lookups_per_seq != 0)
                LOG("max lookups per seq fragment " << max_lookups_per_seq);
        
#if LOOKUP_TABLE
            // determining size of lookup key
            int lookup_key_bits = 1;
            while ((hash_array.size() >> lookup_key_bits) > 5)
            {
                lookup_key_bits += 1;
            }
            hash_lookup_shift = kmer_len * 2 - lookup_key_bits;

            const size_t bucket_count = size_t(1) << lookup_key_bits;
            LOG("creating lookup table with " << bucket_count << " buckets, on average " << (float(hash_array.size()) / bucket_count) << " hashes per bucket");
            hash_lookup_table.resize(bucket_count + 1);

            // figuring out bucket ranges
            size_t hash_idx = 0;
            hash_t last_hash = 0;
            for (size_t bucket_idx = 0; bucket_idx < bucket_count; ++bucket_idx)
            {
                hash_lookup_table[bucket_idx] = hash_idx;
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
            }
            hash_lookup_table[bucket_count] = hash_array.size();
 #endif
        }

        std::pair<tax_t, hash_t> find_hash(hash_t hash, int  default_value ) const
        {
#if LOOKUP_TABLE
            hash_t bucket_idx = hash >> hash_lookup_shift;
            auto first = hash_array.begin() + hash_lookup_table[bucket_idx];
            auto last = hash_array.begin() + hash_lookup_table[bucket_idx + 1];            
#else
            auto first = hash_array.begin();
            auto last = hash_array.end();
#endif
            first = std::lower_bound(first, last, KmerTax(hash, 0));
            if ((first == last) || (hash < first->kmer)){
                return {default_value,default_value};
            }
            else {
                return {first->tax_id, (unique ? first->kmer : default_value)};
            }
            // return { ((first == last) || (hash < first->kmer) ) ? {default_value, default_value}; :  {first->tax_id, first->kmer}; }
        }

        int calculate_lookup_window(int seq_kmers) const
        {
            if ((max_lookups_per_seq == 0) || seq_kmers <= max_lookups_per_seq)
                return 1;

            return (seq_kmers + max_lookups_per_seq - 1) / max_lookups_per_seq;
        }

        Hits operator() (const std::string &seq) const 
        {
            Hits hits;
            int index = 0;
            hash_t min_hash = 0;
            uint64_t min_fnv_hash = 0;

            int seq_kmers = seq.length() - kmer_len + 1;
            const int lookup_window = calculate_lookup_window(seq_kmers);

            Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
                {
                    hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
                    auto fnv_hash = lookup_window == 1 ? 0 : KmerHash::hash_of(hash); 

                    if (index == 0 || fnv_hash < min_fnv_hash)
                    {
                        min_hash = hash;
                        min_fnv_hash = fnv_hash;
                    }

                    index++;

                    if ((index % lookup_window) == 0 || index == seq_kmers) 
                    {
                        auto hit = find_hash(min_hash,0);
                        if (unique && hit.first && hits[hit.first].find(hit.second) == hits[hit.first].end())   {
                            hits[hit.first].emplace(hit.second);
                        }
                        else if (hit.first){
                            hash_t h(hits[hit.first].size() + 1);
                            hits[hit.first].emplace(h);
                        }
                        index = 0;
                    }
                    return true;
                });

            return hits;
        }

    };

    struct KmerTaxIdMatchId
    {
        struct Matches : public std::list<KmerTax>
        {
            operator bool() const { return !empty(); }
        };

    	int seq_id;
        const Matches matches;

    	KmerTaxIdMatchId(int seq_id, const Matches &matches) : seq_id((int)seq_id), matches(matches) {}
    	bool operator < (const KmerTaxIdMatchId &b) const { return seq_id < b.seq_id; }
    };

    struct KmerTaxIdMatcher
    {
        const HashSortedArray &hash_array;
        typedef KmerTaxIdMatchId::Matches Hits;
        int kmer_len;

        KmerTaxIdMatcher(const HashSortedArray &hash_array, int kmer_len) : hash_array(hash_array), kmer_len(kmer_len) {}

        tax_t find_tax_id(hash_t hash) const
        {
            auto first = hash_array.begin();
            auto last = hash_array.end();

            first = std::lower_bound(first, last, KmerTax(hash, 0));
            if ((first == last) || (hash < first->kmer))
                return 0;
            else 
                return first->tax_id;
        }

        Hits operator() (const std::string &seq) const // todo: optimize, move Hits out
        {
            Hits hits;

			Hash<hash_t>::for_all_hashes_do(seq, (int)kmer_len, [&](hash_t hash)
				{
			        hash = seq_transform<hash_t>::min_hash_variant(hash, (int)kmer_len);
                    auto tax_id = find_tax_id(hash);
					if (tax_id != 0)
						hits.push_back(KmerTax(hash, tax_id));

					return true;
				});
            return hits;
        }
    };


    struct TaxMatchId
    {
        int seq_id;
        Hits hits;
        TaxMatchId(int seq_id, const Hits &hits) : seq_id(seq_id), hits(hits)   {}

        bool operator < (const TaxMatchId &b) const { return seq_id < b.seq_id; }
    };


template<class TaxHitsO>
    struct TaxHitsPrinter
    {
        tc::Tax_hits<TaxHitsO> &tax_hits;
        const bool print_counts, compact;
        tc::Spot<TaxHitsO> spot;
        tc::Spot<TaxHitsO> last_spot;

        TaxHitsPrinter(bool print_counts, bool compact, tc::Tax_hits<TaxHitsO> &tax_hits_) : print_counts(print_counts), compact(compact), tax_hits(tax_hits_) {}

        void operator() (const std::vector<Reader::Fragment> &processing_sequences, const std::vector<TaxMatchId> &ids)
        {
            for (auto seq_id : ids) {
                spot.name = processing_sequences[seq_id.seq_id].spotid;
                spot.tax_id.clear();
                if constexpr (TaxHitsO::has_counts()) {
                    spot.counts.clear();
                }
                for_each(spot.name.begin(), spot.name.end(), [](char& c) { if (c == '\n' || c == '\t') c = ' ';});
                for (auto &hit : seq_id.hits) {
                    spot.tax_id.push_back(hit.first);
                    if constexpr (TaxHitsO::has_counts()) {
                        spot.counts.push_back(hit.second.size());
                    }
                }
                if (spot.name == last_spot.name) {
                    last_spot.merge(spot);
                } else {
                    tax_hits.add_row(last_spot);
                    swap(last_spot, spot);
                }
            }
            tax_hits.add_row(last_spot);
        }
    };


    struct TaxPrinter
    {
        IO::Writer &writer;
        const bool print_counts, compact, unique;
        TaxPrinter(bool print_counts, bool compact, IO::Writer &writer, bool unique) : print_counts(print_counts), compact(compact), writer(writer), unique(unique) {}

        void load_uniq_chunk(const std::vector<TaxMatchId> &tm_ids)
        {
            for (auto tm_id : tm_ids){
                for (auto &hit : tm_id.hits){
                    for (auto khit : hit.second){
                        unique_hits[hit.first].emplace(khit);
                    }
                }
            }
        }

        void operator() (const std::vector<Reader::Fragment> &processing_sequences, const std::vector<TaxMatchId> &ids)
        {
            if (unique){
                load_uniq_chunk(ids);
            }
            if (compact)
                print_compact(processing_sequences, ids);
            else
                print(processing_sequences, ids);

            writer.check();
        }

        private:

        void print(const std::vector<Reader::Fragment> &processing_sequences, const std::vector<TaxMatchId> &ids)
        {
            for (auto seq_id : ids)
                print_id_info(processing_sequences, seq_id);
        }

        void print_spotid(const Reader::Fragment &fragment)
        {
            for (auto c : fragment.spotid) 
            {
                if (c == '\t' || c == '\n')
                    c = ' ';
                writer.f() << c;
            }
        }

        void print_id_info(const std::vector<Reader::Fragment> &processing_sequences, const TaxMatchId &seq_id)
        {
            print_spotid(processing_sequences[seq_id.seq_id]);

            for (auto &hit : seq_id.hits) 
            {
                writer.f() << '\t' << hit.first;
                if (print_counts && hit.second.size() > 1)
                    writer.f() << 'x' << hit.second.size();
            }

            writer.f() << std::endl;
        }

        void print_compact(const std::vector<Reader::Fragment> &processing_sequences, const std::vector<TaxMatchId> &ids)
        {
            int from = 0, to = 0;
            for (from = 0; from < ids.size() && processing_sequences.begin()->spotid == processing_sequences[ids[from].seq_id].spotid; from++);
            for (to = (int)ids.size() - 1; to >= 0 && processing_sequences.rbegin()->spotid == processing_sequences[ids[to].seq_id].spotid; to--);

            print_first_and_last(processing_sequences, ids, from, to); // to merge them later correctly

            auto spots = merge_to_spots(processing_sequences, ids, from, to);
            auto reverse_index = spot_tax_reverse_index(spots);

            for (auto &x: reverse_index)
            {
                writer.f() << x.second;
                for (auto &tax : x.first)
                    writer.f() << '\t' << tax;

                writer.f() << std::endl;
            }
        }

        void print_first_and_last(const std::vector<Reader::Fragment> &processing_sequences, const std::vector<TaxMatchId> &ids, int &from, int &to)
        {
            for (int i = 0; i < from; i++)
            {
                writer.f() << '\t';
                print_id_info(processing_sequences, ids[i]);
            }

            for (int i = to + 1; i < ids.size(); i++)
            {
                writer.f() << '\t';
                print_id_info(processing_sequences, ids[i]);
            }
        }

#define INTEGER_SPOTID 0

#if INTEGER_SPOTID
        typedef long long spot_id_t;
#else
        typedef p_string spot_id_t;
#endif
        typedef int tax_id_t;
        typedef std::map<spot_id_t, std::set<tax_id_t> > SpotTaxIds; // todo: unordered map ?

        SpotTaxIds merge_to_spots(const std::vector<Reader::Fragment> &processing_sequences, const std::vector<TaxMatchId> &ids, int &from, int &to)
        {
            SpotTaxIds spots;
            for (int i = from; i <= to; i++)
            {
#if INTEGER_SPOTID
                auto spotid = std::stoll(processing_sequences[ids[i].seq_id].spotid);
#else
                auto spotid = p_string(processing_sequences[ids[i].seq_id].spotid);
#endif
                auto &spot = spots[spotid];

                for (auto &hit : ids[i].hits)
                    spot.insert(hit.first);
            }

            return spots;
        }
        
        typedef std::map< std::set<tax_id_t>, size_t> TaxSpotIds;

        TaxSpotIds spot_tax_reverse_index(const SpotTaxIds &spots)
        {
            TaxSpotIds tax_spots;
            for (auto &spot : spots)
                tax_spots[spot.second]++;

            return tax_spots;
        }

    };

    struct KmerTaxIdPrinter
    {
        IO::Writer &writer;
        int kmer_len = 0;
        KmerTaxIdPrinter(IO::Writer &writer, int kmer_len) : writer(writer), kmer_len(kmer_len){}

    	void operator() (const std::vector<Reader::Fragment> &processing_sequences, const std::vector<KmerTaxIdMatchId> &ids)
    	{
    		for (auto seq : ids)
    		for (auto kmer_tax : seq.matches)
    		    writer.f() << Hash<hash_t>::str_from_hash(kmer_tax.kmer, kmer_len) << '\t' << kmer_tax.tax_id << std::endl;

            writer.check();
    	}
    };


    bool hide_counts = false;
    bool compact = false;


    template<class Options>
    void run_collator(const std::string &filename, IO::Writer &writer, const Config &config)
    {
        //spdlog::stopwatch sw; 

        auto tax_hits = make_unique<tc::Tax_hits<Options>>(true);
        {
            Matcher matcher(hash_array, (int)kmer_len, config.optimization_dbs_max_lookups_per_seq_fragment, config.unique); // todo: move to constructor
            TaxHitsPrinter tc_print(!hide_counts, compact, *tax_hits);

            Job::run_for_matcher(filename, config.spot_filter_file, config.unaligned_only, config.optimization_ultrafast_skip_reader, config.chunk_size,
                [&matcher, &tc_print](const std::vector<Reader::Fragment> &chunk) { 
                    Job::match_and_print<Matcher, TaxHitsPrinter<Options>, TaxMatchId>(chunk, tc_print, matcher);
                    //match_and_print_chunk(chunk, tax_hits, matcher); 
                } );
        }
        // We don't need DB anymore
        hash_array.resize(0); 
        hash_array.shrink_to_fit();
        tax_hits->finalize(); 
        if (config.vectorize) {
            tax_hits->save(filename);
        } else {
            tf::Executor executor(config.num_threads ? config.num_threads : std::thread::hardware_concurrency());
            if (config.compact) {
                auto collated_tax_hits = tax_hits->template collate<tc::tax_hits_options<true, false>>(executor);   
                tax_hits.reset(0);
                collated_tax_hits->group(executor, writer.f());
                collated_tax_hits.reset(0);
            } else {
                auto collated_tax_hits = tax_hits->template collate<Options>(executor);   
                tax_hits.reset(0);
                collated_tax_hits->print(executor, writer.f()); 
                collated_tax_hits.reset(0);
            }
        }
        
    }

    virtual void run(const std::string &filename, IO::Writer &writer, const Config &config) override
    {
        hide_counts = config.hide_counts;
        compact = config.compact;

        if (config.collate || config.vectorize) {
            if (hide_counts)
                run_collator<tc::tax_hits_options<false, false>>(filename, writer, config);
            else if (config.compact)
                run_collator<tc::tax_hits_options<false, false>>(filename, writer, config);
            else
                run_collator<tc::tax_hits_options<false, true>>(filename, writer, config);

        } 
        else 
        {
            if (config.print_kmers_only)
            {
    		    KmerTaxIdMatcher matcher(hash_array, kmer_len);
        		KmerTaxIdPrinter print(writer, kmer_len);
                Job::run_for_matcher(filename, config.spot_filter_file, config.unaligned_only,  config.optimization_ultrafast_skip_reader, config.chunk_size,
                    [&](const std::vector<Reader::Fragment> &chunk) { 
                        Job::match_and_print<KmerTaxIdMatcher, KmerTaxIdPrinter, KmerTaxIdMatchId>(chunk, print, matcher);
                    } );
            }
            else
            {
                Matcher matcher(hash_array, (int)kmer_len, config.optimization_dbs_max_lookups_per_seq_fragment,
                                config.unique); // todo: move to constructor
                TaxPrinter print(!hide_counts, compact, writer, config.unique);

                Job::run_for_matcher(filename, config.spot_filter_file, config.unaligned_only,  config.optimization_ultrafast_skip_reader, config.chunk_size,
                    [&](const std::vector<Reader::Fragment> &chunk) { 
                        Job::match_and_print<Matcher, TaxPrinter, TaxMatchId>(chunk, print, matcher);
                    } );

                if (config.unique){
                    IO::Writer writer_u((!writer.filename.empty()) ? writer.filename + ".uniq" : "");
                    unique_hits.print_uniq_hits(writer_u);
                }
            }
        }
    }
};


DBSJob::Hits DBSJob::unique_hits;


struct DBSBasicJob : public DBSJob
{
    DBSBasicJob(const std::string &dbs)
    {
        kmer_len = DBSIO::load_dbs(dbs, hash_array);
    }
};

#endif

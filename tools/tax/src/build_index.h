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

#include <array>
#include <iostream>
#include "omp_adapter.h"
//#include "kmers.h"
#include "kmer_io.h"
#include "seq_transform.h"
#include "kmer_hash.h"
#include "ready_seq.h"

struct BuildIndex
{
    struct VariableWindowSize
    {
        int suggested_window_size, min_window_size, min_kmers_per_seq;
        VariableWindowSize(int suggested_window_size, int min_window_size, int min_kmers_per_seq) : 
            suggested_window_size(suggested_window_size), min_window_size(min_window_size), min_kmers_per_seq(min_kmers_per_seq){}

        int operator()(int seq_len) const 
        {
            if (seq_len / suggested_window_size < min_kmers_per_seq)
                return std::max(min_window_size, seq_len / min_kmers_per_seq);

            return suggested_window_size;
        }
    };

    template <class Lambda>
    static void process_window(const char *s, int len, int kmer_len, Lambda &&add_kmer)
    {
        if (len < kmer_len)
            return;

        const int THREADS = 32;
        struct ThreadFinding
        {
            KmerHash::hash_of_hash_t min_hash;
            int has_kmer_pos, min_hash_pos;
            ThreadFinding() : min_hash(std::numeric_limits<size_t>::max()), has_kmer_pos(-1), min_hash_pos(-1){}
        };

        std::array<ThreadFinding, THREADS> thread_findings;
        bool has_kmer_found = false;

        #pragma omp parallel num_threads(THREADS)
        for (int i = omp_get_thread_num(); !has_kmer_found && i <= len - kmer_len; i += omp_get_num_threads())
        {
            hash_t kmer = KmerIO::kmer_from(s, i, kmer_len);
            kmer = seq_transform<hash_t>::min_hash_variant(kmer, kmer_len);

            auto thread_id = omp_get_thread_num();

            auto h = KmerHash::hash_of(kmer); // todo: can be optimized

            if (h < thread_findings[thread_id].min_hash)
            {
                thread_findings[thread_id].min_hash = h;
                thread_findings[thread_id].min_hash_pos = i;
            }
        }

        int chosen_kmer_pos = -1;

        {
            auto min_hash = thread_findings[0].min_hash;
            int min_hash_pos = thread_findings[0].min_hash_pos;

            for (int i=1; i<THREADS; i++)
                if (thread_findings[i].min_hash < min_hash)
                {
                    min_hash = thread_findings[i].min_hash;
                    min_hash_pos = thread_findings[i].min_hash_pos;
                }

            chosen_kmer_pos = min_hash_pos;

            if (chosen_kmer_pos < 0)
                throw std::runtime_error("cannot find min hash");
        }

    //  kmers.add_kmer(KmerIO::kmer_from(s, chosen_kmer_pos, kmer_len), tax_id);
        hash_t kmer = KmerIO::kmer_from(s, chosen_kmer_pos, kmer_len);
        kmer = seq_transform<hash_t>::min_hash_variant(kmer, kmer_len);
        add_kmer(kmer, chosen_kmer_pos);
    }

    template <class Lambda>
    static void process_clean_string(p_string p_str, int window_size, int kmer_len, Lambda &&add_kmer)
    {
        for (int start = 0; start <= p_str.len - window_size; start += window_size) // todo: check for integer overflows on very long sequences
        {
            int from = std::max(0, start - (kmer_len - 1));
            int to = std::min(start + window_size, p_str.len);
            process_window(p_str.s + from, to - from, kmer_len, [&](hash_t kmer, int offset) { add_kmer(kmer, offset + from); } );
        }
    }

    template <class Lambda, class WindowSizeLambda>
    static size_t add_kmers(const std::string &filename, WindowSizeLambda &&window_size, int kmer_len, Lambda &&add_kmer)
    {
        Fasta fasta(filename);

        size_t seq_index = 0;
        size_t total_size = 0;
        const int DOT_INTERVAL = 128;

        ReadySeq loading_seq, processing_seq;
        load_sequence(&fasta, &processing_seq);

        while (!processing_seq.seq.empty())
        {
            std::thread loading_thread(load_sequence, &fasta, &loading_seq);

            total_size += processing_seq.seq.size();

            for (auto &clean_string : processing_seq.clean_strings)
                process_clean_string(clean_string, window_size(processing_seq.seq.length()), kmer_len, [&](hash_t kmer, int offset) { add_kmer(kmer); } );

            seq_index++;
            if (seq_index % DOT_INTERVAL == 0)
                std::cerr << ".";

            loading_thread.join();
            swap(processing_seq, loading_seq); // todo: move?
        }

        if (seq_index >= DOT_INTERVAL)
            std::cerr << std::endl;

        return total_size;
    }


    static std::vector<ReadySeq> load_all(const std::string &filename)
    {
        Fasta fasta(filename);

        std::vector<ReadySeq> seqs;
        seqs.reserve(100); // todo: tune

        {
            std::string loaded_seq;

            while (fasta.get_next_sequence(loaded_seq))
            {
                ReadySeq seq;
                seq.seq = loaded_seq;
                seqs.push_back(seq);
            }
        }

        for (auto &seq : seqs)
        {
            SeqCleaner cleaner(seq.seq);
            seq.clean_strings = move(cleaner.clean_strings);
        }

        return seqs;
    }

    static size_t total_len(const SeqCleaner::p_strings &clean_strings)
    {
        size_t len = 0;
        for (auto &s : clean_strings)
            len += s.len;
        return len;
    }

    template <class Lambda>
    static size_t add_kmers_with_markup(const std::string &filename, int window_size, int kmer_len, Lambda &&add_kmer)
    {
        auto seqs = load_all(filename);

        size_t total_size = 0;
        size_t base = 0;
        size_t sum_total_len = 0;
        
        for (auto &processing_seq : seqs)
            sum_total_len += total_len(processing_seq.clean_strings);

        if (sum_total_len > size_t(std::numeric_limits<int>::max()))
            throw std::runtime_error("so long fasta files are not supported in this mode");

        int seq_id = 0;
        for (auto &processing_seq : seqs)
        {
            for (auto &clean_string : processing_seq.clean_strings)
            {
                process_clean_string(clean_string, window_size, kmer_len, [&](hash_t kmer, int offset) { add_kmer(kmer, base + offset, sum_total_len, seq_id); } );
                base += clean_string.len;
            }

            total_size += processing_seq.seq.size();
            seq_id++;
        }

        return total_size;
    }

};

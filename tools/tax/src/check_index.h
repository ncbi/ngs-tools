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

#include <iostream>
#include <chrono>
#include <thread>
#include <array>
#include "omp_adapter.h"
#include "kmer_hash.h"
#include "ready_seq.h"
#include "tax_id_tree.h"
#include "hash.h"
#include "seq_transform.h"
#include <sstream>

template <class Kmers>
struct CheckIndex
{
    static void check_hash(hash_t kmer, int kmer_len, Kmers &kmers, tax_id_t tax_id, std::list<hash_t> &result_hashes)
    {
	    if (kmers.has_kmer_but_not_tax(kmer, tax_id))
		    result_hashes.push_back(kmer);
    }

    static void check_clean_string(Kmers &kmers, p_string p_str, tax_id_t tax_id, int kmer_len)
    {
	    const char *s = p_str.s;
	    int len = p_str.len;

	    if (len < kmer_len)
		    return;

	    const int THREADS = 32;
	    struct ThreadFinding
	    {
		    std::list<hash_t> hashes;
	    };

	    std::array<ThreadFinding, THREADS> thread_findings;

	    #pragma omp parallel num_threads(THREADS)
	    for (int i = omp_get_thread_num(); i <= len - kmer_len; i += omp_get_num_threads())
	    {
		    auto thread_id = omp_get_thread_num();

		    auto kmer = Hash<hash_t>::hash_of(s + i, kmer_len);
	    	kmer = seq_transform<hash_t>::min_hash_variant(kmer, kmer_len);
		    check_hash(kmer, kmer_len, kmers, tax_id, thread_findings[thread_id].hashes);

    #if 0 // ~ x100 times slower!
		    seq_transform<hash_t>::for_all_1_char_variations_do(kmer, kmer_len, [&](hash_t hash)
		    {
			    check_hash(hash, kmer_len, kmers, tax_id, thread_findings[thread_id].hashes);
			    return true;
		    });
    #endif
	    }

	    // todo: multithreaded too ?
	    for (int i=0; i<THREADS; i++)
		    for (auto hash : thread_findings[i].hashes)
			    kmers.add_kmer(hash, tax_id);
    }

    static size_t check_kmers(Kmers &kmers, const std::string &filename, tax_id_t tax_id, int kmer_len) // todo: make generic function
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
			    check_clean_string(kmers, clean_string, tax_id, kmer_len);

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

};


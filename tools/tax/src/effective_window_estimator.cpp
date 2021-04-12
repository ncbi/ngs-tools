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
#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <stdexcept>
#include <random>
#include <iostream>
#include <chrono>
#include "omp_adapter.h"
#include <atomic>

typedef uint64_t hash_t;

#include "log.h"
#include "dbs.h"
#include "ready_seq.h"
//#include "aligns_to_dbss_job.h"
#include "fasta.h"
#include "hash.h"
#include "config_effective_window_estimator.h"
#include "seq_transform.h"
#include "file_list_loader.h"
#include "filename_meta.h"

using DBSAnnot = DBSS::DBSAnnot; 
using DBSAnnotation = DBSS::DBSAnnotation;

using namespace std;
using namespace std::chrono;

const string VERSION = "0.11";

typedef vector<hash_t> HashSortedArray;

bool found(hash_t hash, HashSortedArray &hashes)
{
    return std::binary_search(hashes.begin(), hashes.end(), hash);
}

typedef std::vector<p_string> Windows;

p_string fetch_window(std::vector<ReadySeq> &sequences, int window_len, size_t offset)
{
    size_t shift = 0;
    for (auto &seq : sequences)
    for (auto &s : seq.clean_strings)
    {
        size_t next_shift = shift + size_t(s.len);
        if (next_shift > offset)
        {
            size_t this_string_offset = offset - shift;
            if (this_string_offset + window_len > s.len)
                return p_string(s.s, 0); // out of range

            return p_string(s.s + this_string_offset, window_len);
        }
        shift = next_shift;
    }

    std::cerr << "trying to fetch invalid offset " << offset << " of window len " << window_len << endl;
    throw std::runtime_error("trying to fetch invalid offset");
}

Windows fetch_test_windows(std::vector<ReadySeq> &sequences, int window_len, int windows_to_fetch)
{
    Windows windows;
    auto total_len = ReadySeqOps::total_len(sequences);
//    cerr << "fetch_test_windows total len: " << total_len << endl;
    
    if (total_len < size_t(window_len))
        return windows;
    size_t offset_range = total_len - size_t(window_len);

    std::default_random_engine generator;
    std::uniform_int_distribution<size_t> offset_distribution(0, offset_range);

    for (int i = 0; i < windows_to_fetch; i++)
    {
        auto offset = offset_distribution(generator);
        auto window = fetch_window(sequences, window_len, offset);
        if (window.len > 0)
            windows.push_back(window);
    }

//    cerr << "fetched: " << windows.size() << endl;
    return windows;
}

int identify_windows(const Windows &test_windows, HashSortedArray &hashes, int kmer_len)
{
    std::atomic<int> identified;

    identified = 0;

    const int THREADS = 128;
    #pragma omp parallel num_threads(THREADS) 
    for (int i = omp_get_thread_num(); i < test_windows.size(); i += omp_get_num_threads())
    {
        auto &window = test_windows[i];

        Hash<hash_t>::for_all_hashes_do(window, kmer_len, [&](hash_t hash) {
            hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
            if (found(hash, hashes))
            {
                identified ++;
                return false;
            }
            return true;
        });
    }

    return identified;
}

int estimate_window(const string &filename, int kmer_len, HashSortedArray &hashes, int start_window_len, int max_window_len, int test_windows_count, int min_test_windows_count, float needed_part_identified, float window_upscale_factor, float &best_part_identified)
{
//    const int INFINITE = 1000000000;
    best_part_identified = 0;
    if (hashes.empty())
        return -1;

    auto sequences = load_clean_sequences(filename);
    float window_len = start_window_len;

    while (window_len < max_window_len)
    {
//        cerr << "window len is " << int(window_len) << endl;
        auto test_windows = fetch_test_windows(sequences, int(window_len), test_windows_count);
        if (test_windows.size() <= min_test_windows_count)
            return int(window_len);

        auto identified_windows = identify_windows(test_windows, hashes, kmer_len);
//        cerr << "identified: " << identified_windows << " of " << test_windows.size() << endl;
        best_part_identified = float(identified_windows) / test_windows.size();
        if (best_part_identified >= needed_part_identified)
            return int(window_len);

        window_len *= window_upscale_factor;
    }

    return int(window_len);
}

DBSAnnot find_annotation_by_tax_id(const DBSAnnotation &annotation, int tax_id)
{
    for (auto& annot : annotation) 
        if (annot.tax_id == tax_id)
            return annot;

    return DBSAnnot(tax_id, 0, 0);
}

// todo: move to dbss.h
void load_tax_id_kmers(HashSortedArray &hashes, int tax_id, const string &filename, const DBSAnnotation &annotation)
{
    hashes.clear();
    
    auto annot = find_annotation_by_tax_id(annotation, tax_id);
    if (annot.count == 0)
        return;

    std::ifstream f(filename);
    if (f.fail() || f.eof())
        throw std::runtime_error("cannot open dbss file");

    IO::load_vector_no_size(f, hashes, annot.offset, annot.count);
}

int to_percent(float x)
{
    return int(0.5f + 100*x);
}

int main(int argc, char const *argv[])
{
    Config config(argc, argv);

    auto before = high_resolution_clock::now();

    FileListLoader file_list(config.file_list);

    HashSortedArray hashes;
    int kmer_len = DBSSIO::load_header(config.dbss_filename).kmer_len;

    DBSAnnotation annotation;
    DBSSIO::load_dbs_annotation(config.dbss_filename + ".annotation", annotation);

    for (int i = 0; i < file_list.files.size(); i ++)
    {
        auto file_list_element = file_list.files[i];
        auto tax_id = FilenameMeta::tax_id_from(file_list_element.filename);
        load_tax_id_kmers(hashes, tax_id, config.dbss_filename, annotation);
        float best_part_identified = 0;
        auto window_len = estimate_window(file_list_element.filename, kmer_len, hashes, config.start_window_len, config.max_window_len, config.test_windows_count, config.min_test_windows_count, config.needed_part_identified, config.window_upscale_factor, best_part_identified);
        cout << window_len << '\t' << to_percent(best_part_identified) << '\t' << tax_id << endl;
    }

    cerr << "total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count() << endl;
}

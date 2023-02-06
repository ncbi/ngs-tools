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
#include "dbss.h"
#include "kmer_io.h"
#include "hash.h"
#include "sais/libsais64.h"
#include "bwt/bwt.hpp"

struct DBSSJob : public DBSJob
{

    static char char_from_hash(hash_t hash, int pos)
    {
        hash >>= (2 * pos);
        return Hash<hash_t>::hash_char(hash & 3);
    }

    DBSSJob(const std::string &dbss, const std::string &dbss_tax_list, int num_threads)
    {
        auto dbss_reader = DBSS::make_reader(dbss);
        kmer_len = dbss_reader->header.kmer_len;

        DBSS::DBSAnnotation annotation;
        auto sum_offset = DBSS::load_dbs_annotation(DBSS::DBSAnnot::annotation_filename(dbss), annotation);
        dbss_reader->check_consistency(sum_offset);

        tf::Executor executor(num_threads);

        auto tax_list = DBSS::load_tax_list(dbss_tax_list);
        //DBSS::transform_dbss(dbss_reader, tax_list, annotation);
        DBSS::transform_dbss_matrix(dbss_reader, tax_list, annotation, executor);
        DBSS::read_dbss_matrix(dbss_reader, tax_list, annotation);
        exit(0);

        DBSS::transform_dbss(dbss_reader, tax_list, annotation);
        DBSS::transform_dbss_matrix(dbss_reader, tax_list, annotation, executor);
        //DBSS::transform_dbss_bwt(dbss_reader, tax_list, annotation, 4);
        exit(0);

        DBSS::DBSAnnotation_c annotation_c;
        DBSS::load_dbs_annotation_c("sv.annot", annotation_c);
        //for (const auto& a : annotation_c) {
        //    spdlog::info("{}\t{}\t{}\t{}", a.tax_id, a.offset, a.count, a.is_bm);
        //}
        auto dbss_reader_c = DBSS::make_reader("sv.all");
        HashSortedArray hash_array_c;
        spdlog::info("Loading compressed data");
        DBSS::load_dbss_c(hash_array_c, dbss_reader_c, tax_list, annotation_c, num_threads);
        int c = 0;
        //hash_array_c.resize(100);

        vector<uint8_t> bwt_buffer;
        bwt_buffer.resize(hash_array_c.size());


        vector<int64_t> bwt_tmp_buffer;
        bwt_tmp_buffer.resize(hash_array_c.size() + 1);

        vector<uint64_t> bwt_keys;
        bwt_keys.resize(kmer_len);

        vector<uint64_t> new_hash;
        new_hash.resize(hash_array_c.size());

        vector<uint64_t> test_hash;
        test_hash.resize(hash_array_c.size());

        // BTW Encode
        spdlog::info("Starting encode");
        for (auto kmer_idx = 0; kmer_idx < kmer_len; ++kmer_idx) {
            //spdlog::info("Getting chars for pos {}", kmer_idx);
            for (auto hash_idx = 0; hash_idx < hash_array_c.size(); ++hash_idx) 
                bwt_buffer[hash_idx] = char_from_hash(hash_array_c[hash_idx].kmer, kmer_idx);
            bwt_keys[kmer_idx] = libsais64_bwt(&bwt_buffer[0], &bwt_buffer[0], &bwt_tmp_buffer[0], bwt_buffer.size(), 0, nullptr);

//            auto bwt_key = townsend::algorithm::bwtEncode(bwt_buffer.begin(), bwt_buffer.end());
//            bwt_keys[kmer_idx] = std::distance(bwt_buffer.begin(), bwt_key);

            for (auto hash_idx = 0; hash_idx < hash_array_c.size(); ++hash_idx) 
                new_hash[hash_idx] = Hash<hash_t>::update_hash(bwt_buffer[hash_idx], new_hash[hash_idx]);
        }

        // BTW Decode
        for (auto kmer_idx = 0; kmer_idx < kmer_len; ++kmer_idx) {
            for (auto hash_idx = 0; hash_idx < hash_array_c.size(); ++hash_idx) 
                bwt_buffer[hash_idx] = char_from_hash(new_hash[hash_idx], kmer_idx);
            libsais64_unbwt(&bwt_buffer[0], &bwt_buffer[0], &bwt_tmp_buffer[0], bwt_buffer.size(), nullptr, bwt_keys[(kmer_len - 1) - kmer_idx]);

//            auto decode_key = bwt_buffer.begin();
//            advance(decode_key, bwt_keys[(kmer_len - 1) - kmer_idx]);
//            townsend::algorithm::bwtDecode(bwt_buffer.begin(), bwt_buffer.end(), decode_key);
            for (auto hash_idx = 0; hash_idx < hash_array_c.size(); ++hash_idx) 
                test_hash[hash_idx] = Hash<hash_t>::update_hash(bwt_buffer[hash_idx], test_hash[hash_idx]);
        }

        for (auto hash_idx = 0; hash_idx < hash_array_c.size(); ++hash_idx) {
            string s_old = KmerIO::str_kmer(hash_array_c[hash_idx].kmer, kmer_len);
            string s_new = KmerIO::str_kmer(test_hash[hash_idx], kmer_len);
            if (s_old != s_new)
                throw runtime_error("Round trip failed");

        }
        spdlog::info("{} kmers tested!", hash_array_c.size());
        exit(0);
            //            cout << KmerIO::str_kmer(hash.kmer, kmer_len) << endl;

//        cout << " ------------- new kmers --------------" << endl;
//        for (auto hash : test_hash) 
//            cout << KmerIO::str_kmer(hash, kmer_len) << endl;

/*
        spdlog::info("Loading uncompressed data");
        DBSS::load_dbss(hash_array, dbss_reader, tax_list, annotation, num_threads);

        if (hash_array.size() != hash_array_c.size()) 
            throw runtime_error(fmt::format("Hash size difference! old: {}, new: {}", hash_array.size(), hash_array_c.size()));
        for (size_t i = 0; i < hash_array.size(); ++i) {
            if (hash_array[i].kmer != hash_array_c[i].kmer) 
                throw runtime_error(fmt::format("Element {}: kmer difference", i));
            if (hash_array[i].tax_id != hash_array_c[i].tax_id) 
                throw runtime_error(fmt::format("Element {}: tax_id difference", i));
        }
        spdlog::info("Round-trip test succeeded");
        LOG("Round-trip test passed");
*/

    }
};

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

#include <set>
#include <map>
#include "dbs.h"
#include "log.h"
#include "missing_cpp_features.h"
#include <algorithm>
#include "taskflow/taskflow.hpp"
#include <taskflow/algorithm/sort.hpp>
#include <bm/bm64.h>

#include <bm/bmsparsevec.h>
#include "bm/bmsparsevec_serial.h"
#include "bm/bmserial.h"
#include "sais/libsais64.h"

typedef bm::bvector<> bvector_type;

typedef bm::sparse_vector<uint64_t, bvector_type> sparse_vector_u64_t;

BM_DECLARE_TEMP_BLOCK(TB)
struct DBSS
{
    typedef int tax_id_t;

    struct DBSAnnot
    {
        tax_id_t tax_id;
        size_t count, offset;

        DBSAnnot(int tax_id, size_t count, size_t offset) : tax_id(tax_id), count(count), offset(offset){}

        bool operator < (const DBSAnnot &x) const
        {
            return tax_id < x.tax_id;
        }

        static std::string annotation_filename(const std::string &dbss) { return dbss + ".annotation"; }
    };

    struct DBSAnnot_c
    {
        tax_id_t tax_id;
        size_t count, offset;
        bool is_bm;

        DBSAnnot_c(int tax_id, size_t count, size_t offset, bool is_bm) : tax_id(tax_id), count(count), offset(offset), is_bm(is_bm) {}

        bool operator < (const DBSAnnot_c &x) const
        {
            return tax_id < x.tax_id;
        }

        static std::string annotation_filename(const std::string &dbss) { return dbss + ".annotation"; }
    };


    typedef std::vector<DBSAnnot> DBSAnnotation;
    typedef std::vector<DBSAnnot_c> DBSAnnotation_c;
    typedef std::vector<tax_id_t> TaxList;

    struct DBSSReader
    {
        DBSIO::DBSHeader header;
        const std::string dbss;

        DBSSReader(const std::string &dbss) : dbss(dbss){}
        virtual void check_consistency(size_t sum_offset) = 0;
        virtual void load_kmers(std::vector<hash_t> &hashes, tax_id_t tax_id, const DBSAnnot &annot) = 0;
        virtual void load_kmers_c(std::vector<hash_t> &hashes, tax_id_t tax_id, const DBSAnnot_c &annot, vector<char>& buffer) {}
    };

    struct DBSSFileReader : public DBSSReader
    {
        std::ifstream f;
        DBSSFileReader(const std::string &dbss) : DBSSReader(dbss), f(dbss, std::ios::binary | std::ios::in)
        {
            if (f.fail() || f.eof())
                throw std::runtime_error(std::string("cannot open dbss ") + dbss);

            IO::read(f, header);
        }

        virtual void check_consistency(size_t sum_offset) override
        {
            if (sum_offset != IO::filesize(dbss))
                throw std::runtime_error("inconsistent dbss annotation file");
        }

        virtual void load_kmers(std::vector<hash_t> &hashes, tax_id_t tax_id, const DBSAnnot &annot) override
        {
            hashes.clear();
            IO::load_vector_no_size(f, hashes, annot.offset, annot.count);
        }

        virtual void load_kmers_c(std::vector<hash_t> &hashes, tax_id_t tax_id, const DBSAnnot_c &annot, vector<char>& buffer) override
        {
            hashes.clear();
            if (annot.is_bm) {
        	    f.seekg(annot.offset);
	            if (!f)
		            throw std::runtime_error("load_vector_no_size::cannot seek to offset");
                buffer.resize(annot.count);
                f.read((char*)&buffer[0], annot.count);
                sparse_vector_u64_t sv;
                bm::sparse_vector_deserialize(sv, (const unsigned char*)&buffer[0], TB);
                sv.freeze();
                sparse_vector_u64_t::statistics st;
                sv.calc_stat(&st);
                spdlog::info("Tax_id: {}, kmers: {:L}, memory: {:L}", tax_id, sv.size(), st.memory_used);

                hashes.resize(sv.size());
                sv.decode(&hashes[0], 0, sv.size());

            } else {
                IO::load_vector_no_size(f, hashes, annot.offset, annot.count);
            }
        }
    };

    struct DBSSFolderReader : public DBSSReader
    {
        DBSSFolderReader(const std::string &dbss) : DBSSReader(dbss)
        {
            std::ifstream f(dbss + "/" + "header");
            if (f.fail() || f.eof())
                throw std::runtime_error(std::string("cannot open dbss header ") + dbss);

            f >> header.kmer_len;
        }

        virtual void check_consistency(size_t sum_offset) override {}; // doing nothing at the moment

        virtual void load_kmers(std::vector<hash_t> &hashes, tax_id_t tax_id, const DBSAnnot &annot) override 
        {
            auto filename = tax_id_to_filename(dbss, tax_id);
		    auto kmer_len = DBSIO::load_dbs(filename, hashes);
            if (kmer_len != header.kmer_len)
                throw std::runtime_error(filename + std::string("kmer_len of ") + std::to_string(kmer_len) + " is inconsistent with the header kmer_len of " + std::to_string(header.kmer_len));

            if (annot.count != hashes.size())
                throw std::runtime_error(filename + std::string("kmer number of ") + std::to_string(hashes.size()) + " is inconsistent with the annotation kmer number of " + std::to_string(annot.count));
        };

        static std::string tax_id_to_filename(const std::string &dbss, int tax_id)
        {
            return dbss + "/" + std::to_string(tax_id) + ".db";
        }
    };

    static std::unique_ptr<DBSSReader> make_reader(const std::string &dbss)
    {
        if (IO::file_exists(dbss))
            return make_unique<DBSSFileReader>(dbss);

        if (IO::is_folder(dbss + ".split"))
            return make_unique<DBSSFolderReader>(dbss + ".split");

        throw std::runtime_error(std::string("cannot open dbss ") + dbss);
    }

    static size_t load_dbs_annotation(const std::string &filename, DBSAnnotation &annotation)
    {
        std::ifstream f(filename);
        if (f.fail())
            throw std::runtime_error("cannot open annotation file " + filename);

        size_t offset = sizeof(DBSIO::DBSHeader) + sizeof(size_t);
        tax_id_t prev_tax = 0;

        while (!f.eof())
        {
            DBSAnnot a(0, 0, 0);
            f >> a.tax_id >> a.count;
            //if (f.fail())
            //    break;

            if (!a.count || !a.tax_id)
                continue;
                //throw std::runtime_error("bad annotation format - bad count");

            if (prev_tax >= a.tax_id)
                throw std::runtime_error("bad annotation format - tax less");

            a.offset = offset;
            annotation.push_back(a);
            offset += sizeof(hash_t) *a.count;
            prev_tax = a.tax_id;
        }

        return offset;
    }

    static void load_dbs_annotation_c(const std::string &filename, DBSAnnotation_c &annotation)
    {
        std::ifstream f(filename);
        if (f.fail())
            throw std::runtime_error("cannot open annotation file " + filename);
        size_t file_sz = IO::filesize("sv.all");

        tax_id_t prev_tax = 0;
        bool is_prev_bm = false;
        tax_id_t tax_id;
        size_t count, offset;
        string a_type;
        size_t idx = 0;
        while (!f.eof())
        {
            f >> tax_id >> offset >> count >> a_type;
            //if (f.fail())
            //    break;
            if (idx > 0 && is_prev_bm) {
                auto& a = annotation.back();
                a.count = offset - a.offset;
            }
            if (tax_id == 0)
                break;
            if (!count)
                throw std::runtime_error("bad annotation format - bad count");

            if (prev_tax >= tax_id)
                throw std::runtime_error("bad annotation format - tax less");

            is_prev_bm = a_type == "bm";
            annotation.emplace_back(tax_id, count, offset, is_prev_bm);
            prev_tax = tax_id;
            ++idx;
        }
        //if (idx > 0 && is_prev_bm) {
        //    auto& a = annotation.back();
        //    a.count = file_sz - a.offset;
        //}
        spdlog::info("File size: {}, annot size: {}", file_sz, annotation.size());

    }


    static TaxList load_tax_list(const std::string &filename)
    {
        std::ifstream f(filename);
        if (f.fail())
            throw std::runtime_error("cannot open tax list file");

        TaxList taxes;

        while (!f.eof())
        {
            tax_id_t t = 0;
            f >> t;
            if (f.fail())
                break;

            taxes.push_back(t);
        }

        if (!f.eof())
            throw std::runtime_error("bad tax list file format");

        std::sort(taxes.begin(), taxes.end());
        return taxes;
    }
    template<typename T>
    static uint64_t serialize_vec(ofstream& ofs, T& vec, size_t tax_id)
    {
        vec.optimize(TB);
        vec.freeze();

        sparse_vector_u64_t::statistics st;
        vec.calc_stat(&st);
        spdlog::info("Tax_id: {}, kmers: {:L}, memory: {:L}", tax_id, vec.size(), st.memory_used);

        //bm::sparse_vector_serializer<T> serializer;
        bm::sparse_vector_serial_layout<T> sv_lay;
        bm::sparse_vector_serialize (vec, sv_lay, TB);        
        //serializer.serialize(vec, sv_lay;
        const unsigned char* buf = sv_lay.data();
        ofs.write((char*)&buf[0], sv_lay.size());
        return sv_lay.size();
    }

    static char char_from_hash(hash_t hash, int pos)
    {
        hash >>= (2 * pos);
        return Hash<hash_t>::hash_char(hash & 3);
    }

    static void transform_dbss(std::unique_ptr<DBSSReader> &dbss_reader, const TaxList &tax_list, const DBSAnnotation &annotation)
    {
        LOG("Conversion begin");
        std::vector<hash_t> hashes;
        ofstream ofs("sv.all", ofstream::out | ofstream::binary);
        ofstream ofsa("sv.annot", ofstream::out );
        size_t offset = 0;
        for (auto tax_id : tax_list) {
            for (auto& annot : annotation) {
                if (annot.tax_id == tax_id && annot.count > 0) {
                    dbss_reader->load_kmers(hashes, annot.tax_id, annot);
                    ofsa << annot.tax_id << "\t" << offset << "\t" << hashes.size() << "\t";
                    if (hashes.size() <= 2000) {
                        ofsa << "raw" << "\n";
                        IO::save_vector_data(ofs, hashes);
                        offset += sizeof(uint64_t) * hashes.size();
                    } else {
                        ofsa << "bm" << "\n";        
                        sparse_vector_u64_t sv;
                        sparse_vector_u64_t::back_insert_iterator sv_bi = sv.get_back_inserter();
                        for (auto hash : hashes) 
                            sv_bi.add(hash);
                        sv_bi.flush();
                        offset += serialize_vec(ofs, sv, tax_id);
                    }
                }
            }
        }
        ofsa << 0 << "\t" << offset << "\t\t\n";
        ofs.close();
    }

    static void save_bv(std::ofstream& bfile, int index, bvector_type& bvector) 
    {   
        size_t blob_size = 0;
        try {
            IO::write(bfile, index); 
            if (bvector.empty()) {
                IO::write(bfile, blob_size); 
            } else {
                bvector.optimize(TB);
                typename bvector_type::statistics st1;
                bvector.calc_stat(&st1);
                vector<unsigned char> blob(st1.max_serialize_mem);
                size_t blob_size = bm::serialize(bvector, &blob[0]);
                IO::write(bfile, blob_size); 
                bfile.write((char*)&blob[0], std::streamsize(blob_size));
            }
        } catch (exception& e) {
            spdlog::error("Failed for {}", index);
            throw;
        }
    }

    static void read_bv(std::ifstream& bfile, const size_t len, bvector_type& bvector)
    {
        char* buffer = new char[len];
        bfile.read(buffer, len);
        bm::deserialize(bvector, (unsigned char*)buffer);
        delete [] buffer;
        bvector.freeze();
    }

    struct bit_matrix {
        static const int cSIZE = 65537;
        //static const size_t c48_mask = 0xFFFF000000000000;
        static const size_t c48_mask = 0x0000FFFFFFFFFFFF;
        vector<bvector_type> data;

        bit_matrix(bool init = false) {
            data.resize(cSIZE, bvector_type(bm::BM_GAP));
            if (init)
                for_each(data.begin(), data.end(), [](bvector_type& v) { v.init();});
        }
        void add(size_t kmer) {
            size_t idx = kmer >> 48;
            size_t hash = kmer & c48_mask;
            if (hash == bm::id_max) {
                data[cSIZE - 1].set_bit_no_check(idx);
            } else {
                data[idx].set_bit_no_check(hash);
            }
        }

        bool test(size_t kmer) {
            size_t idx = kmer >> 48;
            size_t hash = kmer & c48_mask;
            if (hash == bm::id_max) 
                return data[cSIZE - 1].test(idx);
            return data[idx].test(hash);               
        }

        mutex m_output_mutex;

        void serialize_all(ofstream& ofs, tf::Executor& executor) 
        {
            for (size_t i = 0; i < cSIZE; ++i)
               save_bv(ofs, i, data[i]); 
            return;

            tf::Taskflow taskflow;
            int step = 1;
            taskflow.for_each_index(0, cSIZE, step, [&ofs, this](int index) {
                size_t blob_size = 0;
                if (data[index].empty()) {
                    const lock_guard<std::mutex> lock(m_output_mutex);
                    IO::write(ofs, index); 
                    IO::write(ofs, blob_size); 
                } else {
                    data[index].optimize();
                    typename bvector_type::statistics st1;
                    data[index].calc_stat(&st1);
                    vector<unsigned char> blob(st1.max_serialize_mem);
                    size_t blob_size = bm::serialize(data[index], &blob[0]);
                    const lock_guard<std::mutex> lock(m_output_mutex);                    
                    IO::write(ofs, index); 
                    IO::write(ofs, blob_size); 
                    ofs.write((char*)&blob[0], std::streamsize(blob_size));
                }
            });
            executor.run(taskflow).wait();
        }

        void serialize(ofstream& ofs, int index) 
        {
            IO::write(ofs, index); 

            size_t blob_size = 0;
            if (data[index].empty()) {
                IO::write(ofs, blob_size); 
            } else {
                data[index].optimize(TB);
                typename bvector_type::statistics st1;
                data[index].calc_stat(&st1);
                vector<unsigned char> blob(st1.max_serialize_mem);
                size_t blob_size = bm::serialize(data[index], &blob[0]);
                IO::write(ofs, blob_size); 
                ofs.write((char*)&blob[0], std::streamsize(blob_size));
            }
        }

        void deserialize(ifstream& ifs) {
            int idx = 0;
            int max_idx = 0;
            size_t sz = 0;
            IO::read(ifs, idx);
            while (!ifs.eof()) {
                max_idx = max<int>(idx, max_idx);
                IO::read(ifs, sz);
                if (sz > 0)
                    read_bv(ifs, sz, data[idx]);
        		ifs.read((char*)&idx, sizeof(idx));
            }
            if (max_idx != cSIZE - 1) {
                spdlog::error("{} != {}", max_idx, cSIZE - 1);
            //    throw runtime_error("Invalid matrix");
            }
        }
        size_t memory_used() {
            typename bvector_type::statistics st;
            size_t mem_used = cSIZE * sizeof(int) + cSIZE * sizeof(size_t);
            for (auto& v : data) {
                if (!v.empty()) {
                    v.calc_stat(&st);
                    mem_used += st.memory_used;
                }
            }
            return mem_used;
        }

    };


    static void test_kmers() {
        static const size_t c48_mask = 0x0000FFFFFFFFFFFF; // First 48 bits
        std::vector<hash_t> test_kmers;
        bm::bvector<> bv;
        bv.init();

        { // read input 
            ifstream ifs("kmers_15829", ofstream::in );
            size_t kmer;
            while (!ifs.eof()) {
                ifs >> kmer;                            
                test_kmers.push_back(kmer);
            }
            spdlog::info("Input done");
            cout << test_kmers.size() << " kmers" << endl;
            for (auto kmer : test_kmers) {
                assert(kmer & c48_mask != bm::id_max);
                bv.set_bit_no_check(kmer & c48_mask);
            }
            spdlog::info("Bit set");

            bv.optimize(TB);
            spdlog::info("Optimized");
            //bv.freeze();
        }
        // check 1
        spdlog::info("Check 1");
        for (auto kmer : test_kmers) {
            if (!bv.test(kmer & c48_mask)) {
                cout << kmer << " failed (bv)" << endl;
            }
        }
        spdlog::info("Saving...");
        bm::SaveBVector("bv_15829", bv);

        spdlog::info("Loading...");
        bm::bvector<> bv1;
        // check 2
        bm::LoadBVector("bv_15829", bv1);
        //bv1.freeze();
        spdlog::info("Check 2");
        for (auto kmer : test_kmers) {
            if (!bv1.test(kmer & c48_mask)) {
                cout << kmer << " failed (bv1)" << endl;
            }
        }

    }

    static void transform_dbss_matrix(std::unique_ptr<DBSSReader> &dbss_reader, const TaxList &tax_list, const DBSAnnotation &annotation,  tf::Executor& executor)
    {
        //tf::Executor executor(num_threads);
        spdlog::info("dbss_matrix_conversion begin");
        test_kmers();
        exit(0);

        bit_matrix m{true};
        std::vector<hash_t> hashes;

        ofstream ofs("bv.matrix", ofstream::out | ofstream::binary);
        size_t offset = 0;

        for (auto tax_id : tax_list) {
            for (auto& annot : annotation) {
                if (annot.tax_id == tax_id && annot.count > 0) {
                    dbss_reader->load_kmers(hashes, annot.tax_id, annot);
/*
                    tf::Taskflow taskflow;
                    int step = 20000;
                    taskflow.for_each_index(0, 65535, step, [&hashes,&step, &m](int i) {
                        static const size_t c48_mask = 0xFFFF000000000000;
                        auto e = i + step;
                        for (auto hash : hashes) {
                            size_t idx = hash;
                            idx >>= 48;
                            if (idx >= i && idx < e) {
                                hash &= (~c48_mask);
                                m[idx].set_bit_no_check(hash);
                            }
                        }
                    });  

                    executor.run(taskflow).wait();
*/          
/*
                    size_t idx = 4455547782370951167 >> 48;
                    spdlog::info("Checking idx {}", idx);
                    for (auto hash : hashes) {
                        if (hash >> 48 != idx)
                            continue;
                        test_hashes.push_back(hash);                            
                        ofs1 << hash << endl;
                        m.add(hash);
                    }
                    m.serialize(ofs, idx);

                    for (auto hash : test_hashes) {
                        if (!m.test(hash)) {
                            spdlog::info("{} not set", hash);
                        }
                    }
                    spdlog::info("m tested");
                    bit_matrix m1;
                    ifstream ifs("bv.matrix", ifstream::in | ifstream::binary);
                    m1.deserialize(ifs);
                    ifs.close();
                    for (auto hash : test_hashes) {
                        if (!m1.test(hash)) {
                            spdlog::info("{} not set", hash);
                        }
                    }
                    spdlog::info("m1 tested");

                    exit(0);
*/                  

                    size_t idx = 0;
                    int step = 1000;
                    while (idx < bit_matrix::cSIZE-1) {
                        auto sz = hashes.size();
                        size_t i = 0;
                        while (i < sz) {
                            size_t curr_idx = hashes[i] >> 48;
                            size_t last_idx = min<int>(bit_matrix::cSIZE-1, idx + step);
                            if (curr_idx >= idx && curr_idx < last_idx) {
                                m.add(hashes[i]);
                                hashes[i] = hashes[sz - 1];
                                --sz;
                                continue;
                            }
                            ++i;
                        }
                        hashes.resize(sz);
                        for (size_t i = 0; i < step; ++i) {
                            if (idx + i >= bit_matrix::cSIZE-1)
                                break;
                            m.serialize(ofs, idx + i);
                            m.data[idx + i].clear(true);
                        }
                        if (idx % 1000 == 0)
                            spdlog::info("saved {} vectors, {}", idx, hashes.size());
                        idx += step;                            
                    }
                    

/*
                    for (size_t idx = 0; idx < bit_matrix::cSIZE-1; ++idx) {
                        size_t last = 0;
                        auto sz = hashes.size();
                        for (size_t i = 0; i < sz; ++i, ++last) {
                            while (i < sz && hashes[i] >> 48 == idx) {
                                m.add(hashes[i]);
                                ++i;
                            }
                            if (i >= sz) break;
                            if (last != i)
                                hashes[last] = hashes[i];
                        }
                        hashes.resize(last);
                        m.serialize(ofs, idx);
                        m.data[idx].clear(true);
                        if (idx % 1000 == 0)
                            spdlog::info("saved {} vectors, {}", idx, hashes.size());
                    }
  */

/*

                    for (size_t idx = 0; idx < bit_matrix::cSIZE-1; ++idx) {
                        size_t last = 0;
                        auto sz = hashes.size();
                        for (size_t i = 0; i < sz; ++i, ++last) {
                            while (i < sz && hashes[i] >> 48 == idx) {
                                m.add(hashes[i]);
                                ++i;
                            }
                            if (i >= sz) break;
                            if (last != i)
                                hashes[last] = hashes[i];
                        }
                        hashes.resize(last);
                        m.serialize(ofs, idx);
                        m.data[idx].clear(true);
                        if (idx % 1000 == 0)
                            spdlog::info("saved {} vectors, {}", idx, hashes.size());
                    }
                    
*/                    
                    m.serialize(ofs, bit_matrix::cSIZE-1);
                }
            }
        }
        spdlog::info("dbss_matrix_conversion end");
        //spdlog::info("dbss_matrix_saving begin");
        //m.serialize_all(ofs, executor);
        //spdlog::info("dbss_matrix_saving end");
        ofs.close();
    }

    static void read_dbss_matrix(std::unique_ptr<DBSSReader> &dbss_reader, const TaxList &tax_list, const DBSAnnotation &annotation)
    {
        spdlog::info("Reading dbss matrix");
        bit_matrix m;
        ifstream ifs("bv.matrix", ifstream::in | ifstream::binary);
        m.deserialize(ifs);
        ifs.close();
        spdlog::info("Reading done");
        
        spdlog::info("Mem used: {:L}", m.memory_used());

        for (auto tax_id : tax_list) {
            for (auto& annot : annotation) {
                if (annot.tax_id == tax_id && annot.count > 0) {
                    std::vector<hash_t> hashes;
                    dbss_reader->load_kmers(hashes, annot.tax_id, annot);
                    for (auto hash : hashes) {
                        if (!m.test(hash)) {
                            spdlog::info("{} not set", hash);
                        }
                    }
                    spdlog::info("{} kmers checked", hashes.size());
                }
            }
        }
    }

    static void transform_dbss_bwt(std::unique_ptr<DBSSReader> &dbss_reader, const TaxList &tax_list, const DBSAnnotation &annotation, int num_threads)
    {
        LOG("BTW Conversion begin");
        std::vector<hash_t> hashes;
        ofstream ofs("sv_bwt.all", ofstream::out | ofstream::binary);
        ofstream ofsa("sv_bwt.annot", ofstream::out );
        size_t offset = 0;
        vector<uint8_t> bwt_buffer;
        vector<int64_t> bwt_tmp_buffer;
        vector<uint64_t> bwt_keys;
        vector<uint64_t> new_hash;
        auto kmer_len = dbss_reader->header.kmer_len;

        for (auto tax_id : tax_list) {
            for (auto& annot : annotation) {
                if (annot.tax_id == tax_id && annot.count > 0) {
                    dbss_reader->load_kmers(hashes, annot.tax_id, annot);
                    //hashes.resize(4);
                    size_t sz = hashes.size();

                    ofsa << annot.tax_id << "\t" << offset << "\t" << sz << "\t";
                    if (hashes.size() <= 2000) {
                        ofsa << "raw" << "\n";
                        IO::save_vector_data(ofs, hashes);
                        offset += sizeof(uint64_t) * sz;
                    } else {
                        ofsa << "bm" << "\n";        
                        bwt_buffer.resize(sz);
                        bwt_tmp_buffer.resize(sz + 1);
                        bwt_keys.resize(kmer_len * 2);
                        new_hash.resize(sz);
/*
                        for (auto kmer_idx = 0; kmer_idx < kmer_len; ++kmer_idx) {
                            for (auto hash_idx = 0; hash_idx < sz; ++hash_idx) 
                                bwt_buffer[hash_idx] = char_from_hash(hashes[hash_idx], kmer_idx);
                            bwt_keys[kmer_idx] = libsais64_bwt_omp(&bwt_buffer[0], &bwt_buffer[0], &bwt_tmp_buffer[0], bwt_buffer.size(), 0, nullptr, num_threads);
                            for (auto hash_idx = 0; hash_idx < sz; ++hash_idx) 
                                new_hash[hash_idx] = Hash<hash_t>::update_hash(bwt_buffer[hash_idx], new_hash[hash_idx]);
                        }
*/                        
                        sparse_vector_u64_t sv;
                        sparse_vector_u64_t::back_insert_iterator sv_bi = sv.get_back_inserter();
                        ofstream of_hashes("kmers1", ofstream::out | ofstream::binary);
                        for (auto hash : hashes) {
                            sv_bi.add(hash);
                            of_hashes << hash << "\n";
                        }
                        sv_bi.flush();
                        of_hashes.close();
                        if ( 1 == 1 ) {
                            sv.optimize(TB);
                            auto num_planes = sv.effective_slices();
                            spdlog::info("Planes: {}, hashes: {:L}", num_planes, hashes.size());
                            for (size_t plane_idx = 0; plane_idx < num_planes; ++plane_idx) {
                                auto plane = sv.slice(plane_idx);
                                if (plane == nullptr)
                                    continue;
                                size_t plane_sz = bwt_buffer.size();
                                fill(bwt_buffer.begin(), bwt_buffer.end(), '0');
                                for (auto en = plane->first(); en.valid(); ++en) {
                                    //plane_sz = *en;
                                    bwt_buffer[*en] = '1';
                                }            

//                                spdlog::info("Plane {}, size {:L}", plane_idx, plane_sz);
                                //spdlog::info("plane {} before", plane_idx);
                                //copy(bwt_buffer.begin(), bwt_buffer.end(), ostream_iterator<char>(std::cout));
                                //cout << endl;
                                bwt_keys[plane_idx] = libsais64_bwt_omp(&bwt_buffer[0], &bwt_buffer[0], &bwt_tmp_buffer[0], bwt_buffer.size(), 0, nullptr, num_threads);

                                //spdlog::info("plane {} after", plane_idx);
                                //copy(bwt_buffer.begin(), bwt_buffer.end(), ostream_iterator<char>(std::cout));
                                //cout << endl;

                                plane->clear(true);
                                for (size_t i = 0; i < plane_sz; ++i)
                                    if (bwt_buffer[i] == '1')
                                        plane->set_bit(i);
                                
                                /*
                                bvector_type new_plane;
                                for (size_t i = 0; i < plane_sz; ++i)
                                    if (bwt_buffer[i] == '1')
                                        new_plane.set_bit(i);

                                if (!plane->equal(new_plane)) 
                                    spdlog::info("Plane {} is not equal", plane_idx);
                                */

                                //auto ins_it = new_plane.inserter();
                                //for (size_t i = 0; i < plane_sz; ++i)
                                //    ins_it = bwt_buffer[i] == '1';

                                //plane->move_from(new_plane);
                                
                            }
                            {
                                sv.optimize(TB);
                                //new_hash.resize(sv.size());
                                auto it = sv.begin();
                                ofstream of_hashes("kmers2", ofstream::out | ofstream::binary);
                                while (it.valid()) {
                                    of_hashes << it.value() << "\n";
                                    it.advance();
                                }
                                of_hashes.close();
/*
                                sv.decode(&new_hash[0], 0, sv.size());
                                ofstream of_hashes("kmers2", ofstream::out | ofstream::binary);
                                for (auto hash : new_hash) {
                                    of_hashes << hash << "\n";
                                }
                                of_hashes.close();
*/                                
                            }

                        }
                        offset += serialize_vec(ofs, sv, tax_id);
                    }
                }
            }
        }
        ofsa << 0 << "\t" << offset << "\t\t\n";
        ofs.close();
        LOG("BTW Conversion done");
    }

    template <class C>
    static void load_dbss(std::vector<C> &hash_array, std::unique_ptr<DBSSReader> &dbss_reader, const TaxList &tax_list, const DBSAnnotation &annotation, int num_threads)
    {
        hash_array.clear();

        size_t total_hashes_count = 0;
        for (auto tax_id : tax_list)
            for (auto& annot : annotation)
                if (annot.tax_id == tax_id)
                    total_hashes_count += annot.count;

        hash_array.reserve(total_hashes_count);
        {        
            std::vector<hash_t> hashes;
            for (auto tax_id : tax_list)
                for (auto& annot : annotation) 
                    if (annot.tax_id == tax_id && annot.count > 0) 
                    {
                        dbss_reader->load_kmers(hashes, tax_id, annot);
                        for (auto hash : hashes) {
                            hash_array.emplace_back(hash, annot.tax_id);
                        }
                    }
        }

        if (hash_array.size() != total_hashes_count)
        {
            std::cerr << hash_array.size() << " of " << total_hashes_count << " loaded " << std::endl;
            throw std::runtime_error("unable to load all kmers");
        }
        
        LOG("dbss parts loaded (" << (total_hashes_count / 1000 / 1000) << "m kmers)");
        if (num_threads > 0) {
            tf::Executor executor(num_threads);
            tf::Taskflow taskflow;
            taskflow.sort(hash_array.begin(), hash_array.end());
            executor.run(taskflow).wait();
        } else {    
            std::sort(hash_array.begin(), hash_array.end()); // todo: parallel sort in new C++
        }
        LOG("dbss parts merged");
    }

    template <class C>
    static void load_dbss_c(std::vector<C> &hash_array, std::unique_ptr<DBSSReader> &dbss_reader, const TaxList &tax_list, const DBSAnnotation_c &annotation, int num_threads)
    {
        hash_array.clear();

        size_t total_hashes_count = 0;
        for (auto tax_id : tax_list)
            for (auto& annot : annotation)
                if (annot.tax_id == tax_id)
                    total_hashes_count += annot.count;

        //hash_array.reserve(total_hashes_count);
        {        
            std::vector<hash_t> hashes;
            vector<char> buffer;
            for (auto tax_id : tax_list)
                for (auto& annot : annotation) 
                    if (annot.tax_id == tax_id && annot.count > 0) 
                    {
                        spdlog::info("Loading tax_id {}", tax_id);
                        dbss_reader->load_kmers_c(hashes, tax_id, annot, buffer);
                        for (auto hash : hashes) {
                            hash_array.emplace_back(hash, annot.tax_id);
                        }
                    }
        }
/*
        if (hash_array.size() != total_hashes_count)
        {
            std::cerr << hash_array.size() << " of " << total_hashes_count << " loaded " << std::endl;
            throw std::runtime_error("unable to load all kmers");
        }
*/        
        LOG("dbss parts loaded_c (" << (total_hashes_count / 1000 / 1000) << "m kmers)");
        if (num_threads > 0) {
            tf::Executor executor(num_threads);
            tf::Taskflow taskflow;
            taskflow.sort(hash_array.begin(), hash_array.end());
            executor.run(taskflow).wait();
        } else {    
            std::sort(hash_array.begin(), hash_array.end()); // todo: parallel sort in new C++
        }
        LOG("dbss parts merged");
    }

};


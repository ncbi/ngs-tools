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

    typedef std::vector<DBSAnnot> DBSAnnotation;
    typedef std::vector<tax_id_t> TaxList;

    struct DBSSReader
    {
        DBSIO::DBSHeader header;
        const std::string dbss;

        DBSSReader(const std::string &dbss) : dbss(dbss){}
        virtual void check_consistency(size_t sum_offset) = 0;
        virtual void load_kmers(std::vector<hash_t> &hashes, tax_id_t tax_id, const DBSAnnot &annot) = 0;
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
            if (f.fail())
                break;

            if (!a.count)
                throw std::runtime_error("bad annotation format - bad count");

            if (prev_tax >= a.tax_id)
                throw std::runtime_error("bad annotation format - tax less");

            a.offset = offset;
            annotation.push_back(a);
            offset += sizeof(hash_t) *a.count;
            prev_tax = a.tax_id;
        }

        return offset;
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
                        for (auto hash : hashes)
                            hash_array.emplace_back(hash, annot.tax_id);
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
};


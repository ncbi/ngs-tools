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

#include <iostream>
#include <chrono>
#include <thread>
#include <array>
#include "omp_adapter.h"
#include "check_index.h"
//#include "kmers.h"
#include "dbs.h"
#include "ready_seq.h"
#include "filename_meta.h"
#include "file_list_loader.h"

#include "ready_seq.h"
#include <sstream>

#include "tax_id_tree.h"
#include "config_identify_tax_ids.h"
//#include "filename_meta.h"

using namespace std;
using namespace std::chrono;


void fail(const std::string &message)
{
    LOG(message);
    throw std::runtime_error(message);
}

struct Kmers
{
    const TaxIdTree &tax_id_tree;

    std::vector<hash_t> storage;
    std::vector<int> tax_ids;

    int kmer_len = 0;

    Kmers(const TaxIdTree &tax_id_tree, const std::string &filename) : tax_id_tree(tax_id_tree)
    {
        kmer_len = DBSIO::load_dbs(filename, storage);
        tax_ids.resize(storage.size());
    }

    void update_kmer(hash_t kmer, tax_id_t tax_id)
    {
        auto pos = find_pos(kmer);
        if (pos == storage.size())
            return;

        auto &at = tax_ids[pos];
        if (at == tax_id) // todo: remove?
            return;

        if (at == 0)
            at = tax_id;
        else
            at = tax_id_tree.consensus_of(tax_id, at);
    }

    size_t find_pos(hash_t kmer) const
    {
        auto it = std::lower_bound(storage.begin(), storage.end(), kmer);
        if (it == storage.end() || *it != kmer)
            return storage.size();

        return std::distance(storage.begin(), it);
    }
};

bool check_if_found_unknown_tax_ids(const TaxIdTree &tax_id_tree, const FileListLoader &file_list)
{
    bool found = false;
    for (auto &file_list_element : file_list.files)
        if (!tax_id_tree.known_id(FilenameMeta::tax_id_from(file_list_element.filename)))
        {
            std::cerr << "unknown tax id: " << FilenameMeta::tax_id_from(file_list_element.filename) << " for file " << file_list_element.filename << endl;
            found = true;
        }

    return found;
}

int main(int argc, char const *argv[])
{
    Config config(argc, argv);
    LOG("identify_tax_ids version 0.12");

    FileListLoader file_list(config.file_list);

    TaxIdTree tax_id_tree;
    TaxIdTreeLoader::load_tax_id_tree(tax_id_tree, config.tax_parents_file);

    if (check_if_found_unknown_tax_ids(tax_id_tree, file_list))
        throw std::runtime_error("unknown tax ids found - exiting");

    Kmers kmers(tax_id_tree, config.db_in_file);
    LOG("kmer len: " << kmers.kmer_len);
    LOG(kmers.storage.size() << " kmers loaded");
    auto before = high_resolution_clock::now();

    size_t total_size = 0;
    for (auto &file_list_element : file_list.files)
    {
        auto tax_id = FilenameMeta::tax_id_from(file_list_element.filename);
        LOG(file_list_element.filesize << "\t" << tax_id << "\t" << file_list_element.filename);
        total_size += CheckIndex<Kmers>::check_kmers(kmers, file_list_element.filename, tax_id, kmers.kmer_len);
        {
            auto seconds_past = std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count();
            if (seconds_past < 1)
                seconds_past = 1;

            size_t megs = total_size/1000000;
            LOG("processed size " << megs << "M = " << (total_size/1000)/seconds_past << "K/sec");
        }
    }

    std::ofstream f(config.out_file);
    IO::save_vector(f, kmers.tax_ids);

    LOG("total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count());
}


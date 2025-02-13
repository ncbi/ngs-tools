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
#include "config_build_index_of_each_file.h"
#include "build_index.h"
#include "file_list_loader.h"
#include "dbs.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.18";

int main(int argc, char const *argv[])
{
    LOG("build_index_of_each_file version " << VERSION);
    Config config(argc, argv);
    LOG("window size: " << config.window_size);
    LOG("kmer len: " << config.kmer_len);

    auto before = high_resolution_clock::now();

    FileListLoader file_list(config.file_list);

#ifdef SINGLETHREADED_BUILD_INDEX
    const int THREADS = 96;

    #pragma omp parallel num_threads(THREADS) 
    for (int i = omp_get_thread_num(); i < file_list.files.size(); i += omp_get_num_threads())
    {
        auto file_list_element = file_list.files[i];

#else
    for (auto &file_list_element : file_list.files)
    {
#endif

        string out_file = file_list_element.filename + config.out_ext;
        if (IO::file_exists(out_file))
        {
            cout << "skipping " << out_file << endl;
            continue;
        }

        std::string summary_file = out_file + ".summary.tsv";

        size_t start = summary_file.rfind("/");
        std::string jf = summary_file.substr(start + 1, summary_file.length() - start);

        size_t end = jf.find(".");
        std::string summary_basename = jf.substr(0, end);

        set<hash_t> kmers;
        std::vector<std::string> summary;
        BuildIndex::add_kmers(file_list_element.filename, BuildIndex::VariableWindowSize(config.window_size, config.min_window_size, config.min_kmers_per_seq), config.kmer_len, summary, [&](hash_t kmer){ kmers.insert(kmer); });

        cout << "writing " << out_file << endl;
        DBSIO::save(out_file, kmers, config.kmer_len);

        cout << "writing summary " << summary_file << endl;
        ofstream fsummary;
        fsummary.open(summary_file);

        for (size_t sit = 0; sit < summary.size(); sit++)
            fsummary << summary_basename << "\t" << summary[sit] << endl;

        fsummary.close();

    }

    LOG("total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count());
}

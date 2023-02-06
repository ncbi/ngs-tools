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
#include <bm/bm64.h>
#include "config_align_to.h"

#include <iostream>
#include <chrono>
#include <thread>
#include <list>
#include "omp_adapter.h"

const std::string VERSION = "0.802";

typedef uint64_t hash_t;

#include "dbs.h"
#include "aligns_to_job.h"
#include "aligns_to_db_job.h" 
#include "aligns_to_dbs_job.h"
#include "aligns_to_dbsm_job.h"
#include "aligns_to_dbss_job.h"
//#include "aligns_to_many_jobs.h"
#include "missing_cpp_features.h"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>
#ifdef __linux__
#include <sys/resource.h>
#endif
#include "taskflow/taskflow.hpp"
#include <taskflow/algorithm/transform.hpp>
#include <chrono>

using namespace std;
using namespace std::chrono;

void test() {

    ifstream ss("text");
    string word; // for storing each word
     // Traverse through all words
    // while loop till we get
    // strings to store in string word
    vector<string> in, out;
    while (ss >> word)
        // print the read word
        in.push_back(word);
    tf::Executor executor(8);
    tf::Taskflow taskflow;
    out.resize(in.size());
    taskflow.transform(in.begin(), in.end(), out.begin(), [](const string& s) {
        std::this_thread::sleep_for(std::chrono::milliseconds(rand() % 50));
        return s;
    });
    executor.run(taskflow).wait();
    for (auto& s : out)
        cout << s << " ";
    exit(0);

}

int main(int argc, char const *argv[])
{
    #ifdef __GLIBCXX__
    std::set_terminate(__gnu_cxx::__verbose_terminate_handler);
    #endif
    //std::locale::global(std::locale("en_US.UTF-8")); // enable comma as thousand separator

    auto stderr_logger = spdlog::stderr_logger_mt("stderr"); // send log to stderr
    spdlog::set_default_logger(stderr_logger);
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v"); // default logging pattern (datetime, error level, error text)

    //test();
    LOG("aligns_to version " << VERSION);
    Config config(argc, argv);
    {
        int num_threads = config.num_threads > 0 ? config.num_threads : omp_get_max_threads();
        LOG("hardware threads: "  << std::thread::hardware_concurrency() << ", omp threads: " << num_threads);
    }
    if (config.num_threads > 0)
        omp_set_num_threads(config.num_threads);

    auto before = high_resolution_clock::now();

    unique_ptr<Job> job;

    if (!config.db.empty())
        job = unique_ptr<DBJob>(new DBJob(config.db));
    else if (!config.dbs.empty())
        job = unique_ptr<DBSBasicJob>(new DBSBasicJob(config.dbs));
    else if (!config.dbsm.empty())
        job = unique_ptr<DBSMJob>(new DBSMJob(config.dbsm));
    else if (!config.dbss.empty())
        job = unique_ptr<DBSSJob>(new DBSSJob(config.dbss, config.dbss_tax_list, config.num_threads));
//    else if (!config.many.empty())
//        job = make_unique<ManyJobs>(config.many);
    else
        Config::fail();

    LOG("loading time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count());
    if (job->db_kmers() > 0)
        LOG("kmers " << job->db_kmers() << " (" << (job->db_kmers() / 1000 / 1000) << "m)");

    for (auto &contig_file : config.contig_files)
    {
        LOG(contig_file);

        IO::Writer writer(config.contig_files.size() == 1 ? config.out : contig_file + config.out);
		try 
		{
	        job->run(contig_file, writer, config);
		}
		catch (std::exception &e)
		{
			LOG(e.what());
			if (config.contig_files.size() == 1)
				throw e;
		}
    }

/*
#ifdef __linux__
    struct rusage r_usage;
    getrusage(RUSAGE_SELF,&r_usage);
    spdlog::info("Total memory: {:L}Kb", r_usage.ru_maxrss);
#endif     
*/

    LOG("total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count());

//    std::exit(0); // dont want to wait for destructors. commented because of some issues with stream flushing
    return 0;
}


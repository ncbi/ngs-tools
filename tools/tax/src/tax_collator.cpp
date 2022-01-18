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

#include "config_tax_collator.h"
#include "tax_collator.hpp"

#include <iostream>
#include <chrono>
#include <thread>
#include <list>

const std::string VERSION = "0.672";
#include <io.h>
#include "missing_cpp_features.h"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>
#ifdef __linux__
#include <sys/resource.h>
#endif


using namespace std;
using namespace std::chrono;

template<class Options>
void run_collator(tf::Executor& executor, bool compact, ostream& os, const string& file) 
{

    //tax_hits_opt::is_compact()
    if (compact) {
        auto tax_hits = std::make_unique<tc::Tax_hits<tc::tax_hits_options<true, false>>>(true);
        tax_hits->init(file);
        auto collated_tax_hits = tax_hits->template collate<tc::tax_hits_options<true, false>>(executor);   
        tax_hits.reset(0);
        collated_tax_hits->group(executor, os);
        collated_tax_hits.reset(0);
        /*
        auto tax_hits = std::make_unique<tc::Tax_hits<Options>>(executor, true);
        tax_hits->init(file);
        auto collated_tax_hits = tax_hits->template collate<tc::tax_hits_options<true, false>>();   
        tax_hits.reset(0);
        //collated_tax_hits->save(file + ".ccc");        
        collated_tax_hits->group(os);
        collated_tax_hits.reset(0);
        */
    } else {
        auto tax_hits = std::make_unique<tc::Tax_hits<Options>>(true);
        tax_hits->init(file);
        auto collated_tax_hits = tax_hits->template collate<Options>(executor);   
        tax_hits.reset(0);
        collated_tax_hits->print(executor, os); 
    }
}


int main(int argc, char const *argv[])
{
    #ifdef __GLIBCXX__
    std::set_terminate(__gnu_cxx::__verbose_terminate_handler);
    #endif

    std::locale::global(std::locale("en_US.UTF-8")); // enable comma as thousand separator

    auto stderr_logger = spdlog::stderr_logger_mt("stderr"); // send log to stderr
    spdlog::set_default_logger(stderr_logger);
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v"); // default logging pattern (datetime, error level, error text)
    
    LOG("tax_collator " << VERSION);
    LOG("hardware threads: "  << std::thread::hardware_concurrency());
    Config config(argc, argv);
    tf::Executor executor{(size_t)config.num_threads};

    auto before = high_resolution_clock::now();

    for (auto &contig_file : config.contig_files)
    {
        LOG(contig_file);

        IO::Writer writer(config.contig_files.size() == 1 ? config.out : contig_file + config.out);
        try 
        {
            
            if (config.hide_counts)
                run_collator<tc::tax_hits_options<false, false>>(executor, config.compact, writer.f(), contig_file) ;
            else if (config.compact)
                run_collator<tc::tax_hits_options<false, false>>(executor, config.compact, writer.f(), contig_file);
            else
                run_collator<tc::tax_hits_options<false, true>>(executor, config.compact, writer.f(), contig_file);
            
        }
        catch (std::exception &e)
        {
            LOG(e.what());
            if (config.contig_files.size() == 1)
                throw e;
        }
    }

#ifdef __linux__
    struct rusage r_usage;
    getrusage(RUSAGE_SELF,&r_usage);
    spdlog::info("Total memory: {:L}Kb", r_usage.ru_maxrss);
#endif     

    LOG("total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count());

//    std::exit(0); // dont want to wait for destructors. commented because of some issues with stream flushing
    return 0;
}


#include "config_align_to.h"

#include "kmer_map.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <list>

using namespace std;
using namespace std::chrono;

#ifdef _OPENMP
   #include <omp.h>
#else
   inline int omp_get_max_threads() { return 0; }
#endif

#define COUNT_MATCHES 0

const string VERSION = "0.35";

#include "merge.h"

typedef uint64_t hash_t;

#include "dbs.h"
#include "aligns_to_job.h"
//#include "database_io.h"
#include "aligns_to_db_job.h"
#include "aligns_to_dbs_job.h"
#include "aligns_to_dbss_job.h"

void print_current_time()
{
	auto t = std::time(nullptr);
	auto timeinfo = std::localtime(&t);
	const int BUFFER_SIZE = 256;
	char buffer[BUFFER_SIZE];
	strftime(buffer, BUFFER_SIZE, "%m/%d/%Y %H:%M:%S", timeinfo);
	cerr << "time is " << buffer << endl;
}

int main(int argc, char const *argv[])
{
    std::set_terminate(__gnu_cxx::__verbose_terminate_handler);
    
    cerr << "aligns_to version " << VERSION << endl;
    cerr << "hardware threads: "  << std::thread::hardware_concurrency() << ", omp therads: " << omp_get_max_threads() << std::endl;
    Config config(argc, argv);

    print_current_time();
    auto before = high_resolution_clock::now();

    Job *job = nullptr;

    if (!config.db.empty())
        job = new DBJob(config);
    else if (!config.dbs.empty())
        job = new DBSBasicJob(config);
    else if (!config.dbss.empty())
        job = new DBSSJob(config);
    else
        Config::fail();

    cerr << "loading time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;
    if (job->db_kmers() > 0)
        cerr << "kmers " << job->db_kmers() << " (" << (job->db_kmers() / 1000 / 1000) << "m)" << endl;

    if (!config.contig_files.empty())
	{
        for (auto &filename : config.contig_files)
            {
                cerr << filename << endl;
                before = high_resolution_clock::now();
                {
                    ofstream out_f(filename + ".matches");
                    out_f.flush(); // ?
                    job->run(filename, out_f);
                }

                auto processing_time = std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count();
                cerr << "processing time (sec) " << processing_time << endl;
            }
    }
    else
	{
        if (config.contig_file.empty())
            throw std::runtime_error("contig file(s) is empty");

        cerr << config.contig_file << endl;
        job->run(config.contig_file, cout);
    }

    cerr << "total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;

    exit(0); // dont want to wait for destructors
    return 0;
}


#include <iostream>
#include <chrono>
#include <set>
#include <map>
#include "config_build_freq_db_conv.h"
#include "../seq_loader.h"
#include "../file_list_loader.h"
#include "../hash.h"
#include <omp.h>

using namespace std;
using namespace std::chrono;

typedef int hash_t;
typedef std::map<hash_t, size_t> Frequences;

const string VERSION = "0.10";

void save_frequences(const string &filename, const Frequences &freqs)
{
    std::ofstream f(filename);
    f.flush();

    f << freqs.size() << endl;
    for (auto &it : freqs)
        f << it.first << " " << it.second << endl;

    if (f.fail())
        throw std::runtime_error("failed to save frequences");
}

string get_out_filename(const string &filename, int kmer_len)
{
    return filename + ".conv." + std::to_string(kmer_len) + ".mer";
}

hash_t get_kmer_from(const char *s, int from, int len)
{
	return Hash<hash_t>::hash_of(s + from, len);
}

void build_freq_db(const string &filename, int kmer_len)
{
	Frequences freqs;

    SeqLoader::for_every_clean_sequence_do(filename, [&](p_string line)
    {
		for (int kmer_from = 0; kmer_from <= line.len - kmer_len; kmer_from++)
		{
            auto kmer = get_kmer_from(line.s, kmer_from, kmer_len);
			freqs[kmer]++;
		}
    });

    string out_file = get_out_filename(filename, kmer_len);
    save_frequences(out_file, freqs);
};

void build_freq_db_list(const Config &config)
{
	FileListLoader file_list(config.file_list);
    int file_number = 0;

    const int THREADS = 16;

	#pragma omp parallel num_threads(THREADS)
    for (int file_number = omp_get_thread_num(); file_number < int(file_list.files.size()); file_number += THREADS)
    {
        auto &file = file_list.files[file_number];
        cerr << file_number << " of " << file_list.files.size() << " loading file " << file.filename << endl;
        build_freq_db(file.filename, config.kmer_len);
    }
}

int main(int argc, char const *argv[])
{
    try
    {
		cerr << "build_freq_db_conv version " << VERSION << endl;
		Config config(argc, argv);
		cerr << "kmer len: " << config.kmer_len << endl;

		auto before = high_resolution_clock::now();

        build_freq_db_list(config);

		cerr << "total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;
        return 0;
    }
    catch ( exception & x )
    {
        cerr << x.what() << endl;
		cerr << "exit 3" << endl;
		return 3;
    }
    catch ( string & x )
    {
        cerr << x << endl;
		cerr << "exit 4" << endl;
		return 4;
    }
    catch ( const char * x )
    {
        cerr << x << endl;
		cerr << "exit 5" << endl;
		return 5;
    }
    catch ( ... )
    {
        cerr << "unknown exception" << endl;
//		cerr << "exit 6" << endl;
		return 6;
    }
}

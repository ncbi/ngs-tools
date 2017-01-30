#include <iostream>
#include <chrono>

//#include <thread>
#include <set>
//#include <list>
//#include <iomanip>
//#include <algorithm>
//#include <map>
//#include <omp.h>
//#include <cmath>
#include "config_build_freq_db.h"
#include "seq_loader.h"
#include "freq_db.h"
#include "freq_db_io.h"
#include "file_list_loader.h"
#include "hash.h"
#include <omp.h>
//#include "io.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.12";

string get_out_filename(const string &filename, int kmer_len)
{
    return filename + ".freq" + std::to_string(kmer_len) + ".amino";
}

hash_t get_kmer_from(const char *s, int from, int len)
{
//	return Hash<hash_t>::hash_of(s + from, len);
	return Hash<hash_t>::hash_of_amino(s + from, len);
//	return seq_transform<hash_t>::min_hash_variant2(kmer, len);
}

void get_window_from_to(int &from, int &to, int kmer_from, int line_len, int kmer_len, int window_r)
{
	from = std::max(0, kmer_from - window_r);
	to = std::min(line_len - kmer_len + 1, kmer_from + window_r + 1);
}

template <class Lambda>
void do_for_left_and_right_sides(p_string line, int kmer_from, int kmer_len, int window_r, Lambda &&lambda)
{
	int from, to;
	get_window_from_to(from, to, kmer_from, line.len, kmer_len, window_r);
	for (int pos = from; pos <= kmer_from - kmer_len; pos++)
        lambda(pos);
    
	for (int pos = kmer_from + kmer_len; pos < to; pos++)
        lambda(pos);
}

bool same_frame(int pos1, int pos2)
{
    const int FRAME_NUCLEOTIDES = 3;
    return ( (pos1 - pos2) % FRAME_NUCLEOTIDES ) == 0;
}

void update_dictionary_window(p_string line, hash_t kmer, int kmer_from, int kmer_len, int window_r, PredictedFrequences &pred_freqs)
{
	set<hash_t> added_kmers; // todo: think. remove ?
	auto &pred_freqs_of_kmer = pred_freqs[kmer];
//	cout << "update window kmer_from: " << kmer_from << " from " << from << " to " << to << endl;

    do_for_left_and_right_sides(line, kmer_from, kmer_len, window_r, [&](int pos)
		{
			if (pos + kmer_len > line.len) // todo: remove
				throw "pos + kmer_len > int(line.size())";

#if 0
            if (!same_frame(pos, kmer_from)) // todo: remove or optimoze
                return;
#endif

            auto with_kmer = get_kmer_from(line.s, pos, kmer_len);

			if (added_kmers.find(with_kmer) == added_kmers.end())
            {
    			added_kmers.insert(with_kmer);
	    		pred_freqs_of_kmer[with_kmer] ++;
    //			cout << kmer << " " << with_kmer << endl;
            }
		});
}

void build_freq_db(const string &filename, int kmer_len, int window_r)
{
	Frequences freqs;
	PredictedFrequences pred_freqs;
//    const int kmer_len = config.kmer_len;

    SeqLoader::for_every_clean_sequence_do(filename, [&](p_string line)
    {
		hash_t prev_kmer = 0; // todo: ?
		for (int kmer_from = 0; kmer_from <= line.len - kmer_len; kmer_from++)
		{
            auto kmer = get_kmer_from(line.s, kmer_from, kmer_len);
			if (prev_kmer == kmer) // todo: check for prev prev too ? ACACACAC ?
				continue;

			prev_kmer = kmer;

			freqs.add(kmer);
			update_dictionary_window(line, kmer, kmer_from, kmer_len, window_r, pred_freqs);
		}
    });

#if 0
	for (auto &it : pred_freqs)
		Freq::to_freq(it.second, freqs.freqs[it.first]);
#else
	for (auto &it : pred_freqs)
		Freq::to_freq(it.second, Freq::calculate_total(it.second));
#endif

	Freq::to_freq(freqs.freqs, freqs.total);

    string out_file = get_out_filename(filename, kmer_len);
//    cerr << "saving to " << out_file << endl;
    FreqDBIO::save_frequences(out_file, freqs, pred_freqs, kmer_len, window_r);
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
        build_freq_db(file.filename, config.kmer_len, config.window_r);
    }
}

int main(int argc, char const *argv[])
{
    try
    {
		cerr << "build_freq_db version " << VERSION << endl;
		Config config(argc, argv);
		cerr << "kmer len: " << config.kmer_len << endl;
		cerr << "window_r : " << config.window_r << endl;

		auto before = high_resolution_clock::now();

//        if (!config.fasta_filename.empty())
  //          build_freq_db(config.fasta_filename, config.kmer_len, config.window_r);
    //    else
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

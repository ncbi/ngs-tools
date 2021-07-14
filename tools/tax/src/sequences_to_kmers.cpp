#include <iostream>
#include <chrono>
#include <atomic>
#include <set>
#include <map>
#include <random>
#include "config_sequences_to_kmers.h"
#include "file_list_loader.h"
#include "seq_transform.h"
#include "io.h"

typedef uint64_t hash_t;

#include "hash.h"
#include "fasta.h"

#include "omp_adapter.h"

using namespace std;
using namespace std::chrono;

#define DO_PARALLEL_PER_FILE 1
#define SAVE_IN_CURRENT_FOLDER 0

typedef vector<hash_t> Kmers;

void add_kmers(Kmers &kmers, const string &seq, int kmer_len)
{
    Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
    {
        hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
        kmers.push_back(hash);
        return true;
    });
}

void save(const string &filename, const Kmers &kmers)
{
//    cout << "saving to " << filename << endl;
    std::ofstream f(filename, std::ios::out | std::ios::binary);
    IO::write(f, kmers.size()); // todo: save_vector
    for (auto &kmer : kmers)
        IO::write(f, kmer);
}

string nodir(const string &filename)
{
    auto pos = filename.find_last_of('/');
    if (pos == string::npos)
        return filename;

    return filename.substr(pos + 1);
}

string save_file(const string &filename)
{
#if SAVE_IN_CURRENT_FOLDER
    return string("./") + nodir(filename) + ".profile";
#else
    return string("./") + filename + ".profile";
#endif
}

void to_kmers(const string &filename, int kmer_len)
{
    Fasta fasta(filename);
    string seq;
    Kmers kmers;
    kmers.reserve(10000); // todo: tune

    while (fasta.get_next_sequence(seq))
        add_kmers(kmers, seq, kmer_len);

    save(save_file(filename), kmers);
}

void to_kmers(const Config &config)
{
	FileListLoader file_list(config.file_list);

    atomic<int> counter;
    counter = 0;
#if DO_PARALLEL_PER_FILE
   	#pragma omp parallel for // num_threads(96)
#endif
    for (int file_number = 0; file_number < int(file_list.files.size()); file_number ++)
    {
        auto &file = file_list.files[file_number];
        to_kmers(file.filename, config.kmer_len);
        cerr << '.';
        counter++;
        if ((counter % 1000) == 0)
            cout << endl << counter << endl;
    }
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	auto before = high_resolution_clock::now();

    try
    {
        to_kmers(config);
    }
    catch (std::exception &ex)
    {
        cerr << ex.what() << endl;
    }

	cerr << "total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;
}

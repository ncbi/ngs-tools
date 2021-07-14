#include <iostream>
#include <chrono>
#include <atomic>
#include <set>
#include <map>
#include <random>
#include "io.h"
#include "config_print_profile.h"

typedef uint64_t hash_t;

#include "hash.h"

using namespace std;

vector<hash_t> load_profile(const string &filename)
{
    std::ifstream f(filename, std::ios::in | std::ios::binary);
    if (!f.good())
        throw std::runtime_error(string("cannot load profile ") + filename);

    vector<hash_t> kmers;
    IO::load_vector(f, kmers);
    return kmers;
}

void print_profile(const string &filename)
{
    auto kmers = load_profile(filename);
    for (auto kmer : kmers)
        cout << Hash<hash_t>::str_from_hash(kmer, 32) << endl;
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

    print_profile(config.filename);
}

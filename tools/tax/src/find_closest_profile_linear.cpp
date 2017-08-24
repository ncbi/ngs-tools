#include <iostream>
#include <chrono>
#include <set>
#include <map>
#include "config_find_closest_profile_linear.h"
#include "file_list_loader.h"
#include "io.h"
#include <algorithm>

typedef uint64_t hash_t;

using namespace std;
using namespace std::chrono;

struct Profile
{
    string filename;
    vector<hash_t> kmers;
};

typedef std::vector<Profile> Profiles;

void load_profile(const string &filename, Profile &profile)
{
    std::ifstream f(filename, std::ios::in | std::ios::binary);
    if (!f.good())
        throw std::runtime_error(string("cannot load profile ") + filename);

    profile.filename = filename;
    IO::load_vector(f, profile.kmers);
}

void load_profiles(const string &file_list_name, Profiles &profiles)
{
	FileListLoader file_list(file_list_name);
    profiles.resize(file_list.files.size());
    cout << "loading " << file_list.files.size() << " profiles" << endl;
    for (int file_number = 0; file_number < int(file_list.files.size()); file_number ++)
    {
        auto &file = file_list.files[file_number];
        load_profile(file.filename, profiles[file_number]);
//        load_profile(file.filename + ".profile", profiles[file_number]);
    }
    cout << "loaded" << endl;
}

struct ProfileSimilarity
{
    double operator() ( const Profile &a, const Profile &b) const
    {
        if (a.kmers.size() != b.kmers.size())
            throw std::runtime_error("ProfileSimilarity:: a.kmers.size() != b.kmers.size()");

        if (a.kmers.empty())
            throw std::runtime_error("ProfileSimilarity:: a.kmers.empty()");

        double sum = 0;
        for (int i=0; i < a.kmers.size(); i++)
            if (a.kmers[i] == b.kmers[i])
                sum++;

        return sum / a.kmers.size();
    }
};

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	auto before = high_resolution_clock::now();

    Profiles profiles;
    load_profiles(config.file_list, profiles);

    Profile profile;
    load_profile(config.profile_file, profile);

    ProfileSimilarity sim;

    struct ComparisonResult
    {
        double sim = 0;
        string filename;

        ComparisonResult() = default;
        ComparisonResult(double sim, const string &filename) : sim(sim), filename(filename){}

        bool operator < (const ComparisonResult &x) const
        {
            return sim > x.sim;
        }
    };

    typedef std::vector<ComparisonResult> ComparisonResults;

    ComparisonResults res(profiles.size());

    const int THREADS = 32;
   	#pragma omp parallel for num_threads(THREADS)
    for (int i = 0; i < res.size(); i++)
        res[i] = ComparisonResult(sim(profile, profiles[i]), profiles[i].filename);

    std::sort(res.begin(), res.end());

    for (int i = 0; i < config.top_count && i < res.size(); i++)
        cout << res[i].sim << " " << res[i].filename << endl;

	cerr << "total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;
}

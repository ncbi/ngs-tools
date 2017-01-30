#include "similarity_of.h"
#include "config_similarity_of.h"
#include "file_list_loader.h"

void compress_list(const Config &config)
{
    struct Result
    {
        string filename;
        double sim;
        Result(const string &filename = "", double sim = 0) : filename(filename), sim(sim){}

        bool operator < (const Result &r) const
        {
            return sim > r.sim;
        }
    };

	FileListLoader file_list(config.freq_file_list);
    int file_number = 0;
    vector<Result> results(file_list.files.size());

    const int THREADS = 16;

    int a_kmer_len = 0, a_window_r = 0;
	Frequences a_freqs;
	PredictedFrequences a_pred_freqs;
    FreqDBIO::load_frequences(config.a_freq_filename, a_freqs, a_pred_freqs, a_kmer_len, a_window_r);

	#pragma omp parallel num_threads(THREADS)
    for (int file_number = omp_get_thread_num(); file_number < int(file_list.files.size()); file_number += THREADS)
    {
        auto &file = file_list.files[file_number];
        cerr << file_number << " of " << file_list.files.size() << " loading file " << file.filename << endl;
//        auto sim = similarity_of(config.a_freq_filename, file.filename + ".freq9.amino");
        auto sim = similarity_of(a_pred_freqs, a_freqs.freqs, file.filename + ".freq9.amino");
        results[file_number] = Result(file.filename, sim);
        {
            auto &r = results[file_number];
            cerr << r.sim << '\t' << r.filename << endl;
        }
    }

    std::sort(results.begin(), results.end());

    cout << "results:" << endl;
    for (auto &r : results)
        cout << r.sim << '\t' << r.filename << endl;
}

int main(int argc, char const *argv[])
{
    try
    {
		Config config(argc, argv);

		auto before = high_resolution_clock::now();

        if (!config.b_freq_filename.empty())
            cout << similarity_of(config.a_freq_filename, config.b_freq_filename) << endl;
        else
            compress_list(config);

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

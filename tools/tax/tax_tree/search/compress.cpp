#include "compress.h"
#include "config_compress.h"

void compress_list(const Config &config)
{
    struct Result
    {
        string filename;
        double bits;
        double intra_part;
        Result(const string &filename = "", double bits = 0, double intra_part = 0) : filename(filename), bits(bits), intra_part(intra_part) {}

        bool operator < (const Result &r) const
        {
//            return intra_part < r.intra_part; // 
            return bits < r.bits;
        }
    };

	FileListLoader file_list(config.freq_file_list);
    int file_number = 0;
    vector<Result> results(file_list.files.size());

//    for (auto &file : file_list.files)
//    for (int file_number = 0; file_number < int(file_list.files.size()); file_number++)
    const int THREADS = 12;

	#pragma omp parallel num_threads(THREADS)
    for (int file_number = omp_get_thread_num(); file_number < int(file_list.files.size()); file_number += THREADS)
    {
        auto &file = file_list.files[file_number];
        cerr << file_number << " of " << file_list.files.size() << " loading file " << file.filename << endl;
        auto bits = compress(file.filename + ".freq9.amino", config.filename);
        results[file_number] = Result(file.filename, bits.bits, bits.intra_part);
        {
            auto &r = results[file_number];
            cerr << r.bits << '\t' << r.intra_part << '\t' << r.filename << endl;
        }
    }

    cout << results.size() << endl;
    for (auto &r : results)
        cout << r.bits << endl;

    std::sort(results.begin(), results.end());

    cout << "results:" << endl;
    for (auto &r : results)
        cout << r.bits << '\t' << r.intra_part << '\t' << r.filename << endl;
}

int main(int argc, char const *argv[])
{
    try
    {
		Config config(argc, argv);

		auto before = high_resolution_clock::now();

        if (!config.freq_filename.empty())
        {
            auto r = compress(config.freq_filename, config.filename);
            cout << "bits: " << r.bits << " " << r.intra_part << endl;
        }
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

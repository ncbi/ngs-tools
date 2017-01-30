#include "config_split_to_sequences.h"
#include "fasta.h"
#include "file_list_loader.h"
#include <omp.h>

using namespace std;
//using namespace std::chrono;

void save_sequence(const string &filename, const string &desc, const string &seq)
{
    ofstream f(filename);
    f.flush();
    f << ">" << desc << endl;
    f << seq;
}

bool is_refseq(const string &desc)
{
    return desc.find("|ref|") != string::npos;
}

void split(const string &fasta_filename)
{
    Fasta fasta(fasta_filename);
    string seq;
    const int MAX_SEQUENCES_TO_SAVE = 1000;

    for(int i=0; i < MAX_SEQUENCES_TO_SAVE && fasta.get_next_sequence(seq); )
    {
        string new_filename = fasta_filename + "." + std::to_string(i);
        string desc = fasta.sequence_description();
        if (!is_refseq(desc))
            continue;

        save_sequence(new_filename, desc, seq);

        #pragma omp critical
        {
            cout << Fasta::filesize(new_filename) << "\t" << new_filename << endl;
        }

        i++;
    }
}

void split_list(const Config &config)
{
	FileListLoader file_list(config.file_list);
    int file_number = 0;

    const int THREADS = 1;

	#pragma omp parallel num_threads(THREADS)
    for (int file_number = omp_get_thread_num(); file_number < int(file_list.files.size()); file_number += THREADS)
    {
        auto &file = file_list.files[file_number];
        cerr << file_number << " of " << file_list.files.size() << " loading file " << file.filename << endl;
        split(file.filename);
    }
}

int main(int argc, char const *argv[])
{
    try
    {
		Config config(argc, argv);

//		auto before = high_resolution_clock::now();

        split_list(config);

	//	cerr << "total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;
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

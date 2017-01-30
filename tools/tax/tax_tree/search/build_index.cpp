#include "build_index.h"

const string VERSION = "0.35";

size_t weight(size_t kmers_count)
{
	return kmers_count*16;
}

bool is_eukariota(const string &filename)
{
	return filename.find("/Eukaryota/") != string::npos;
}

bool is_virus(const string &filename)
{
	return filename.find("/Viruses/") != string::npos;
}

int main(int argc, char const *argv[])
{
    try
    {
		cerr << "build_index version " << VERSION << endl;
		ConfigBuildIndex config(argc, argv);
		cerr << "window divider: " << config.window_divider << endl;
		cerr << "kmer len: " << config.kmer_len << endl;

		print_current_time();
		auto before = high_resolution_clock::now();

		FileListLoader file_list(config.file_list);

		TaxIdTree tax_id_tree;
		TaxIdTreeLoader::load_tax_id_tree(tax_id_tree, config.tax_parents_file);

		Kmers kmers(tax_id_tree);
		size_t total_size = 0;
		for (auto &file_list_element : file_list.files)
		{
			auto window_size = calculate_window_size(file_list_element.filesize, is_eukariota(file_list_element.filename), is_virus(file_list_element.filename), config.window_divider);
			auto tax_id = tax_id_from(file_list_element.filename);
			cerr << file_list_element.filesize << "\t" << window_size << "\t" << tax_id << "\t" << file_list_element.filename << endl;
			total_size += add_kmers(kmers, file_list_element.filename, tax_id, window_size, config.kmer_len);
			{
				auto seconds_past = std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count();
				if (seconds_past < 1)
					seconds_past = 1;

				size_t megs = total_size/1000000;
				cerr << "processed size " << megs << "M = " << (total_size/1000)/seconds_past << "K/sec, kmers: " << kmers.storage.size()/1000 << "K, compression rate " << total_size/std::max(size_t(1), weight(kmers.storage.size())) << endl;
			}
		}

		print_kmers(kmers, config.kmer_len);

		cerr << "total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count() << endl;
//		cerr << "total time (msec) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;

		exit(0); // dont want to wait for KmerMaps destructors
        return 0;
    }
    catch ( exception & x )
    {
        cerr << x.what() << endl;
//		cerr << "exit 3" << endl;
		return 3;
    }
    catch ( string & x )
    {
        cerr << x << endl;
//		cerr << "exit 4" << endl;
		return 4;
    }
    catch ( const char * x )
    {
        cerr << x << endl;
//		cerr << "exit 5" << endl;
		return 5;
    }
    catch ( ... )
    {
        cerr << "unknown exception" << endl;
//		cerr << "exit 6" << endl;
		return 6;
    }
}

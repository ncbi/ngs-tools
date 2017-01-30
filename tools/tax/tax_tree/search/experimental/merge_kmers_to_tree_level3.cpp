#include <iostream>
#include <chrono>
#include <set>
#include <map>
#include "config_merge_kmers_to_tree_level3.h"
#include "../seq_loader.h"
#include "../file_list_loader.h"

#include "features.h"
#include "feature_tree.h"
#include "feature_tree_io.h"
#include "feature_tree_builder.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.10";

void merge_kmers_to_tree(TreeBuilder<SimilarityOf3Mers> &builder, const string &filename, int data_frame, int window_len, int between_windows_step, const Tree &tree_3mers, const Tree &tree_9mers)
{
    SeqLoader::for_every_clean_sequence_do(filename, [&](p_string line)
    {
        int count = 0;
		for (int window_from = data_frame; window_from <= line.len - window_len; window_from += between_windows_step, count++)
        {
            builder.insert(calculate_features_for_9mers(line.s + window_from, window_len, tree_3mers, tree_9mers));
            if (count % 10 == 0)
                cerr << ".";
        }
    });
};

void merge_kmers_to_tree_list(const Config &config)
{
    Tree tree_3mers; 
    TreeIO::load_tree(tree_3mers, config.tree_3mer_file);
    SimilarityOf3Mers sim3(tree_3mers);

    Tree tree_9mers; 
    TreeIO::load_tree(tree_9mers, config.tree_9mer_file);
    SimilarityOfNMers simn(tree_9mers);

    TreeBuilder<SimilarityOfNMers> builder(simn);

    if (TreeIO::file_exists(config.tree_file))
        TreeIO::load_tree(builder.tree, config.tree_file);

	FileListLoader file_list(config.file_list);
    int file_number = 0;

    for (int file_number = 0; file_number < int(file_list.files.size()); file_number ++)
    {
        auto &file = file_list.files[file_number];
        cerr << file_number << " of " << file_list.files.size() << " loading file " << file.filename << endl;
        merge_kmers_to_tree(builder, file.filename, config.data_frame, config.window_len, config.between_windows_step, tree_3mers, tree_9mers);

        cout << "---------------------- tree nodes:" << builder.tree.nodes.size() << endl;
        cerr << "---------------------- tree nodes:" << builder.tree.nodes.size() << endl;

        cout << "saving: levels: " << config.levels_to_save << endl;
//        TreeIO::print_tree(builder.tree, config.levels_to_save);
        TreeIO::save_tree(builder.tree, config.tree_file, config.levels_to_save);
        TreeIO::load_tree(builder.tree, config.tree_file);
        cout << "loaded: " << endl;
//        TreeIO::print_tree(builder.tree);

        cout << "loaded  -------------- tree nodes:" << builder.tree.nodes.size() << endl;
        cerr << "loaded  -------------- tree nodes:" << builder.tree.nodes.size() << endl;
    }

    TreeIO::print_tree(builder.tree);
}

int main(int argc, char const *argv[])
{
    try
    {
		cerr << "merge_kmers_to_tree_level3 version " << VERSION << endl;
		Config config(argc, argv);

		auto before = high_resolution_clock::now();

        merge_kmers_to_tree_list(config);

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

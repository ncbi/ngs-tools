#include <iostream>
#include <chrono>
#include <set>
#include <map>
#include "config_compare_1_level.h"
#include "../seq_loader.h"
#include "../file_list_loader.h"

#include "feature_tree_builder.h"
#include "feature_tree_io.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.10";

ClosestNode choose_closest(const vector<ClosestNode> &nodes)
{
    if (nodes.empty())
        throw std::runtime_error("choose_closest:: node.empty()");

    ClosestNode closest = nodes[0];
    for (int i = 1; i<int(nodes.size()); i++)
        if (nodes[i].similarity > closest.similarity)
            closest = nodes[i];

    return closest;
}

Features calculate_features(const string &filename, const Tree &tree_3mers)
{
    Features result;
    SimilarityOf1Mers sim_1_mers;

    SeqLoader::for_every_clean_sequence_do(filename, [&](p_string line)
    {
        const int KMER_LEN = 3;
        const int FRAME_LEN = 3;

        vector<ClosestNode> closest_nodes(FRAME_LEN);
		for (int window_from = 0; window_from <= line.len - KMER_LEN - FRAME_LEN; window_from += FRAME_LEN)
        {
            for (int frame = 0; frame < FRAME_LEN; frame++) // todo: parallel for
            {
                auto fs = FeatureOperations::calculate_features(line.s + window_from + frame, KMER_LEN);
                closest_nodes[frame] = find_closest_final_node(tree_3mers, fs, sim_1_mers);
            }

            auto closest = choose_closest(closest_nodes);
            if (!tree_3mers.valid_node(closest.id))
                std::cerr << "!valid_node(closest.id)";
            else
                result[closest.id] += closest.similarity;
        }
    });

    result.norm = FeatureOperations::norm(result); // todo: type normalized ?
    return result;
};

void compress_1_level(const Config &config)
{
    Tree tree_3mers; 
    TreeIO::load_tree(tree_3mers, config.tree_3mer_file);

    auto fa = calculate_features(config.file_a, tree_3mers);
    auto fb = calculate_features(config.file_b, tree_3mers);

    SimilarityOfNMers simn(tree_3mers);
    auto levels = simn.calculate_leveled_similarity(fa, fb);
    cout << "levels: " << endl;
    double sum = 0, impact = 0.5;
    for (auto &l : levels)
    {
        cout << l << endl;
        sum += l * impact;
        impact *= 0.5;
    }

    cout << "total sum : " << sum << endl;
    cout << "verification " << simn(fa, fb) << endl;
}

int main(int argc, char const *argv[])
{
    try
    {
		cerr << "compress_1_level version " << VERSION << endl;
		Config config(argc, argv);

		auto before = high_resolution_clock::now();

        compress_1_level(config);

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

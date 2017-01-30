#ifndef FEATURE_TREE_BUILDER_H_INCLUDED
#define FEATURE_TREE_BUILDER_H_INCLUDED

#include "feature_tree.h"
#include <set>
#include <omp.h>

static bool almost_identical(double similarity)
{
    return similarity >= 0.99;
}

#if 0
static Features extend_features_to_parents(const Features &features, const Tree &tree)
{
    Features target; // = features;
    for (auto &n : features)
    {
        auto count = n.second;
        tree.for_node_and_parents_do( n.first, [&](const Tree::Node &node)
        {
            target[node.id] += count;
            count *= 0.5; // todo: think
        });
    }

    target.norm = FeatureOperations::norm(target);
    return target;
}
#endif

struct ClosestNode
{
    Tree::node_id id;
    double similarity;

    ClosestNode(Tree::node_id id = Tree::INVALID_NODE, double similarity = 0) : id(id), similarity(similarity){}
};

struct SimilarityOf1Mers
{
    double operator ()(const Features &a, const Features &b) const
    {
        return FeatureOperations::get_similarity(a, b);
    }
};

struct SimilarityOf3Mers // todo: rename to SimilarityOfNMers
{
    const Tree &tree_3mers;
    SimilarityOf3Mers(const Tree &tree_3mers) : tree_3mers(tree_3mers){}

    double operator ()(const Features &a, const Features &b) const
    {
        auto sims = calculate_leveled_similarity(a, b);
        double sum_similarity = 0;
        double sim_impact = 0.5;
        for (auto sim : sims)
        {
            sum_similarity += sim * sim_impact;
            sim_impact *= 0.5;
        }

        return sum_similarity;
    }

    std::vector<double> calculate_leveled_similarity(const Features &a, const Features &b) const
    {
        auto a_levels = get_levels(a, tree_3mers);
        auto b_levels = get_levels(b, tree_3mers);
        int max_level = std::max(get_max_level(a_levels), get_max_level(b_levels));
        if (max_level == 0)
            throw std::runtime_error("SimilarityOf3Mers:: max_level is 0");

        std::vector<double> sims;
        sims.reserve(max_level);
        for (int level = 1; level <= max_level; level++) // todo: parallel for ?
        {
            //std::cout << "on level " << level << std::endl; //<< " similarity is " << sim << std::endl;
            //std::cout << "a_levels[level]" << std::endl;
            //print_features(a_levels[level]);
            //std::cout << "b_levels[level]" << std::endl;
            //print_features(b_levels[level]);

            auto sim = FeatureOperations::get_similarity(a_levels[level], b_levels[level]);
  //          sum_similarity += sim;
//            sum_similarity += sim * sim_impact;
            sims.push_back(sim);
        }

        return sims;
    //    return sum_similarity/max_level;
    }

    typedef std::map<int, Features> Levels;

    static int get_max_level(const Levels &ls)
    {
        int max_level = 0;
        for (auto &l : ls)
            max_level = std::max(max_level, l.first);

        return max_level;
    }
    
    static Levels get_levels(const Features &fs, const Tree &tree)
    {
        Levels levels;

        for (auto &f : fs)
        {
            auto node_id = f.first;
            auto count = f.second;

            int level = tree.get_height_level(node_id);
//            std::cout << "for node " << node_id << " level: " << level << std::endl;

            tree.for_node_and_parents_do(node_id, [&](const Tree::Node &node)
            {
                if (level > 0) // we don't need 0, dont use it anyway
                {
                    auto &level_features = levels[level];
                    level_features[node.id] += count;
                }
                if (level < 0)
                    throw std::runtime_error("get_levels:: level < 0");

                level--;
            });

            for (auto &l : levels)
                l.second.norm = FeatureOperations::norm(l.second);
        }

        return levels;
    }
};

typedef SimilarityOf3Mers SimilarityOfNMers;

template <class SimilarityOf>
static ClosestNode find_closest_final_node(const Tree &tree, const Features &features, SimilarityOf &&similarity_of)
{
    if (tree.empty())
        return ClosestNode(Tree::INVALID_NODE, 0);

    auto node_id = tree.head_id;
    auto similarity = similarity_of(tree.get_node(node_id).features, features); //FeatureOperations::get_similarity(tree.get_node(node_id).features, features);

    while (true)
    {
//        std::cout << node_id << std::endl;
        const auto &node = tree.get_node(node_id);

        if (node.final())
        {
  //          std::cout << "o";
            return ClosestNode(node.id, similarity);
        }

        auto a_sim = similarity_of(tree.get_node(node.a).features, features);//FeatureOperations::get_similarity(tree.get_node(node.a).features, features);
        auto b_sim = similarity_of(tree.get_node(node.b).features, features);//FeatureOperations::get_similarity(tree.get_node(node.b).features, features);

        if (a_sim > b_sim )
        {
            node_id = node.a;
            similarity = a_sim;
        }
        else 
        {
            node_id = node.b;
            similarity = b_sim;
        }
    }
}

template <class SimilarityOf>
struct TreeBuilder
{
    Tree tree;
    const SimilarityOf &similarity_of;
    TreeBuilder(const SimilarityOf &similarity_of) : similarity_of(similarity_of){}

    void insert(const Features &features)
    {
        if (features.empty())
            return;

//        auto closest = find_closest_node(features);
//        auto closest = find_closest_node_flexible(features, std::min(size_t(1000000000), tree.nodes.size()));
        auto closest = find_closest_final_node(tree, features, similarity_of);

        if (!tree.valid_node(closest.id))
        {
            tree.set_head(features);
            return;
        }

//        auto &closest_node = tree.get_node(closest.id);
//        if (tree.get_node(closest.id).final() && almost_identical(closest.similarity))
//            tree.merge_nodes(closest.id, features);
//        else
            tree.insert_neighbour(closest.id, features);
    }
};


template <class Lambda>
static void for_every_kmer_of(const char *s, int kmer_len, int window_len, int step, Lambda &&f)
{
//    #pragma omp parallel for num_threads(3) // todo: think
    for (int i = 0; i <= window_len - kmer_len; i += step)
        f(s + i);
}

static Features calculate_features_for_3mers(const char *_s, int window_len, const Tree &tree_3mers)
{   
    Features result;
//    SimilarityOf3Mers sim(tree_3mers);
    SimilarityOf1Mers sim_1_mers;

    const int KMER_LEN = 3;
    const int FRAME_LEN = 3;

    for_every_kmer_of(_s, KMER_LEN, window_len, FRAME_LEN, [&](const char *kmer)
    {
//        const double MIN_SIMILARITY = 0.5; // todo: think

        auto fs = FeatureOperations::calculate_features(kmer, KMER_LEN);
        auto closest = find_closest_final_node(tree_3mers, fs, sim_1_mers);
//        if (closest.similarity > MIN_SIMILARITY)
        {
            if (!tree_3mers.valid_node(closest.id))
                throw std::runtime_error("calculate_features_for_3mers:: !tree_3mers.valid_node(closest.id)");

//            #pragma omp critical
            {
                result[closest.id] += closest.similarity;
            }
        }
    });

    result.norm = FeatureOperations::norm(result); // todo: type normalized ?
    return result;
}

#define FAST_FINAL_NODE_SEARCH 1

static Features calculate_features_for_9mers(const char *_s, int window_len, const Tree &tree_3mers, const Tree &tree_9mers)
{   
    Features result;
    const int INTERNAL_WINDOW_LEN = 9; // todo: remove hardcoded values
    const int FRAME_LEN = 3;
    const int INTERNAL_KMER_LEN = 3;

#if FAST_FINAL_NODE_SEARCH
    SimilarityOf1Mers sim_n_mers; //(tree_3mers);
#else
//    SimilarityOfNMers sim_n_mers(tree_3mers);
#endif

    for_every_kmer_of(_s, INTERNAL_WINDOW_LEN, window_len, FRAME_LEN, [&](const char *window_start)
    {
        auto fs_3mers = calculate_features_for_3mers(window_start, INTERNAL_WINDOW_LEN, tree_3mers);
        auto closest = find_closest_final_node(tree_9mers, fs_3mers, sim_n_mers);

        {
//                throw std::runtime_error("calculate_features_for_9mers:: !tree_9mers.valid_node(closest.id)");
            if (!tree_9mers.valid_node(closest.id))
                std::cerr << "!tree_9mers.valid_node(closest.id)";
            else
            {
//                #pragma omp critical
                {
                    result[closest.id] += closest.similarity;
                }
            }
        }

    });

    result.norm = FeatureOperations::norm(result); // todo: type normalized ?
    return result;
}


static Features calculate_features_for_27mers(const char *_s, int window_len, const Tree &tree_3mers, const Tree &tree_9mers, const Tree &tree_27mers)
{   
    Features result;
    const int INTERNAL_WINDOW_LEN = 27; // todo: remove hardcoded values
    const int FRAME_LEN = 3;
    const int INTERNAL_KMER_LEN = 3;

#if FAST_FINAL_NODE_SEARCH
    SimilarityOf1Mers sim_9_mers; //(tree_9mers);
#else
//    SimilarityOfNMers sim_9_mers(tree_9mers); // todo: swith to 1 mers?
#endif

    for_every_kmer_of(_s, INTERNAL_WINDOW_LEN, window_len, FRAME_LEN, [&](const char *window_start)
    {
        auto fs_9mers = calculate_features_for_9mers(window_start, INTERNAL_WINDOW_LEN, tree_3mers, tree_9mers);
        auto closest = find_closest_final_node(tree_27mers, fs_9mers, sim_9_mers);

        {
//                throw std::runtime_error("calculate_features_for_9mers:: !tree_9mers.valid_node(closest.id)");
            if (!tree_27mers.valid_node(closest.id))
                std::cerr << "!tree_27mers.valid_node(closest.id)";
            else
            {
//                #pragma omp critical
                {
                    result[closest.id] += closest.similarity;
                }
            }
        }

    });

    result.norm = FeatureOperations::norm(result); // todo: type normalized ?
    return result;
}

#endif
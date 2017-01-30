#include <iostream>
#include <chrono>
#include <thread>
#include <sstream>
#include <set>
#include <omp.h>
#include <algorithm>
#include "compress.h"
#include "file_tree.h"
#include "config_compress_tree.h"

using namespace std;
using namespace std::chrono;

const string EXTENSION = ".freq9.amino";
const int INF_BITS = 1000000;

struct Choice
{
    const Tree::Node *node;
    string full_name, full_folder;
    double bits;
    int level;
    Choice(const Tree::Node *node, const string &full_name, const string &full_folder) : node(node), full_name(full_name), full_folder(full_folder), bits(INF_BITS), level(split(full_folder, '/').size()){}

    bool operator < (const Choice &x) const
    {
//        return bits < x.bits;
        return level * 100 + bits < x.level * 100 + x.bits;
    }
};

typedef vector<Choice> Choices;

string get_full_name(const Tree::Node &node, const string &parent)
{
    string full_name = node.name; //parent + "/" + node.name;
    if (!parent.empty())
        full_name = parent + "/" + full_name;

    return full_name;
}

Choice get_merged_choice(const Tree::Node &node, const string &parent_folder)
{
    string full_name = get_full_name(node, parent_folder);
        
    if (final_node(node))
        return Choice(&node, full_name, parent_folder);

    if (node.subnodes.size() > 1)
        return Choice(&node, full_name + "/merged" + EXTENSION, full_name);

    const Tree::Node &first = node.subnodes.begin()->second;
    return get_merged_choice(first, full_name);
}

Choices get_choices(const Choice &parent)
{
    Choices cs;
    for (auto &subnode: parent.node->subnodes)
        cs.push_back(get_merged_choice(subnode.second, parent.full_folder));

    return cs;
}

void find_bits_for(Choices &cs, const string &fasta_file)
{
    for (auto &c : cs)
    {
        c.bits = compress(c.full_name, fasta_file);
        cout << c.bits << " - bits for " << "\t" << c.full_name << endl;
    }
}

Choice find_best_choice(const Choices &choices)
{
    Choice best_choice(nullptr, "", "");
    for (auto &c : choices)
        if (c.bits < best_choice.bits)
            best_choice = c;

    return best_choice;
 }

void compress_tree(Tree::Node &node, const string &fasta_file)
{
//    Choice best_choice(&node, "", "");
    set<Choice> rating, rating_not_tried_deeper; 
    rating_not_tried_deeper.insert(Choice(&node, "", ""));

    while (true)
    {
        if (rating_not_tried_deeper.empty())
            break;

        Choice best_choice = *rating_not_tried_deeper.begin();
        rating_not_tried_deeper.erase(best_choice);

        cout << "trying deeper for " << best_choice.bits << "\t" << best_choice.full_name << endl << endl;

        auto choices = get_choices(best_choice);
        find_bits_for(choices, fasta_file);
        for (auto &c : choices)
        {
            rating.insert(c);
            rating_not_tried_deeper.insert(c);
        }

        auto best_of = find_best_choice(choices);
//        if (!choices.empty() && best_of.bits > rating.begin()->bits * 1.3)
//        {
//            cout << "stopping " << best_of.bits << " " << best_of.full_name << endl;
//            break;
////            return best_choice.bits;
//        }

//        best_choice = best_of;
//        cout << "best is " << best_choice.bits << "\t" << best_choice.full_name << endl << endl;
    }

    cout << "rating" << endl;
    for (auto &c : rating)
        cout << c.bits << "\t" << c.full_name << endl;
}

double compress_tree_simple(Tree::Node &node, const string &fasta_file)
{
    Choice best_choice(&node, "", "");

    while (true)
    {
        auto choices = get_choices(best_choice);
        find_bits_for(choices, fasta_file);
        auto best_of = find_best_choice(choices);
        if (best_of.bits > best_choice.bits)
        {
            cout << best_choice.bits << " " << best_choice.full_name << endl;
            return best_choice.bits;
        }

        best_choice = best_of;
        cout << "best is " << best_choice.bits << "\t" << best_choice.full_name << endl << endl;
    }
}

int main(int argc, char const *argv[])
{
    try
    {
		Config config(argc, argv);

		auto before = high_resolution_clock::now();
        auto files = load_files_list(config.files_list, EXTENSION);

        Tree tree;
        build_tree(tree, files);
//        print_tree(tree.root);
        cout << "tree loaded" << endl;
        compress_tree(tree.root, config.fasta_file);

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

#include <iostream>
#include <chrono>
#include <thread>
#include <sstream>
#include <set>
#include <omp.h>
#include <algorithm>
#include "freq_db.h"
#include "freq_db_io.h"
#include "file_list_loader.h"

#include "config_merge_freq_db.h"

using namespace std;
using namespace std::chrono;

typedef std::list<string> Strings;

struct Tree
{
    struct Node
    {
        string name;
        map<string, Node> subnodes;

        Node(const string &name = "") : name(name){}
    } root;

    void insert(const Strings &folders)
    {
        if (folders.empty())
            throw "insert empty folders";

        if (root.name.empty())
            root.name = *folders.begin();
        else
            if (root.name != *folders.begin())
                throw "root.name != *folders.begin()";

        auto folder_it = folders.begin();
        folder_it++;
        insert_node(folder_it, folders, root);
    }

    void insert_node(Strings::const_iterator folder_it, const Strings &folders, Node &node)
    {
        if (folder_it == folders.end())
            return;

        if (node.subnodes.find(*folder_it) == node.subnodes.end())
            node.subnodes[*folder_it] = Node(*folder_it);

        auto &to = node.subnodes[*folder_it];
        folder_it++;
        insert_node(folder_it, folders, to);
    }

//    Tree() : root("."){}
};

void split(const string &s, char delim, list<string> &elems) 
{
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim))
        elems.push_back(item);
}

list<string> split(const string &s, char delim) 
{
    list<string> elems;
    split(s, delim, elems);
    return elems;
}

void build_tree(Tree &tree, const Strings &files)
{
    for (auto &name : files)
    {
        auto folders = split(name, '/');
        tree.insert(folders);
    }
}

void print_tree(Tree::Node &node, const string &parent = "")
{
    string full_name = node.name; //parent + "/" + node.name;
    if (!parent.empty())
        full_name = parent + "/" + full_name;

    cout << full_name << endl;
    for (auto &child : node.subnodes)
        print_tree(child.second, full_name);
}

void save_merge_file(const string &merge_file, const Strings &to_merge)
{
    if (to_merge.size() < 2)
        throw "to_merge.size() < 2";

    cout << "saving to " << merge_file << endl;
    for (auto &s : to_merge)
        cout << "->\t" << s << endl;

    Frequences freqs;
    PredictedFrequences pred_freqs;
    int kmer_len = 0, window_r = 0;
    auto to_merge_it = to_merge.begin();

    FreqDBIO::load_frequences(*to_merge_it, freqs, pred_freqs, kmer_len, window_r);
    to_merge_it++;
    for (; to_merge_it != to_merge.end(); to_merge_it++)
    {
        Frequences freqs2;
        PredictedFrequences pred_freqs2;
        int kmer_len2 = 0, window_r2 = 0;
        FreqDBIO::load_frequences(*to_merge_it, freqs2, pred_freqs2, kmer_len2, window_r2);
        if (kmer_len != kmer_len2 || window_r != window_r2)
            throw "kmer_len != kmer_len2 || window_r != window_r2";

        Freq::add(freqs.freqs, freqs2.freqs);
        Freq::add(pred_freqs, pred_freqs2);
    }

    Freq::normalize(freqs.freqs, to_merge.size());
    Freq::normalize(pred_freqs, to_merge.size());

    FreqDBIO::save_frequences(merge_file, freqs, pred_freqs, kmer_len, window_r);
}

const string EXTENSION = ".freq9.amino";

bool final_node(Tree::Node &node)
{
    return node.subnodes.empty();
}

string merge_tree(Tree::Node &node, const string &parent = "")
{
    string full_name = node.name; //parent + "/" + node.name;
    if (!parent.empty())
        full_name = parent + "/" + full_name;

    if (final_node(node))
        return full_name;

    Strings files_to_merge;
//    auto files_to_merge = merge_tree(node.subnodes);
    for (auto &subnode: node.subnodes)
        files_to_merge.push_back(merge_tree(subnode.second, full_name));

    if (files_to_merge.empty())
        throw "files_to_merge.empty()!";

    if (files_to_merge.size() == 1)
        return *files_to_merge.begin();

    full_name += "/merged" + EXTENSION;
    save_merge_file(full_name, files_to_merge);
    return full_name;
}

bool file_exists(const string &filename)
{
    ifstream f(filename);
    return f.good();
}

Strings load_files_list(const string &files_list, const string &ext)
{
    Strings names;
	FileListLoader file_list(files_list);
    for (auto &f : file_list.files)
    {
        string filename = f.filename + ext;
        if (file_exists(filename))
            names.push_back(filename);
        else
            cerr << "file skipped: " << filename << endl;
    }

    return names;
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
        print_tree(tree.root);
        merge_tree(tree.root);

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

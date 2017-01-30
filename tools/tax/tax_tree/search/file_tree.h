#ifndef FILE_TREE_INCLUDED
#define FILE_TREE_INCLUDED

#include <list>
#include <string>
#include "file_list_loader.h"

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

static void split(const string &s, char delim, list<string> &elems) 
{
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim))
        elems.push_back(item);
}

static std::list<string> split(const string &s, char delim) 
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

bool final_node(const Tree::Node &node)
{
    return node.subnodes.empty();
}

//bool file_exists(const string &filename)
//{
//    ifstream f(filename);
//    return f.good();
//}

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

#endif
#ifndef FEATURE_TREE_IO_H_INCLUDED
#define FEATURE_TREE_IO_H_INCLUDED

#include "feature_tree.h"
#include <fstream>

struct TreeIO
{
    template <class IDCast = LowerLevelFeatureIndex>
    static void print_tree(const Tree &tree, int levels = 1000000000)
    {
        return print_tree<IDCast>(tree, tree.head_id, levels, 0);
    }

private:
    template <class IDCast = LowerLevelFeatureIndex>
    static void print_tree(const Tree &tree, Tree::node_id id, int levels, int indent)
    {
        const auto &node = tree.get_node(id);
        if (indent >= levels)
            return;

        for (int i=0; i < indent; i++)
            std::cout << '\t';

        std::cout << "id " << id << '\t' << "norm " << int(node.features.norm + 0.5) << '\t';
        print_features<IDCast>(node.features);

        if (!node.final())
        {
            if (tree.get_node(node.a).parent != node.id)
                throw std::runtime_error("tree.get_node(a).parent != node.id");
            if (tree.get_node(node.b).parent != node.id)
                throw std::runtime_error("tree.get_node(b).parent != node.id");

            print_tree<IDCast>(tree, node.a, levels, indent + 1);
            print_tree<IDCast>(tree, node.b, levels, indent + 1);
        }
    }

    template <class X>
    static void write(std::ofstream &f, X x)
    {
		f.write((char*)&x, sizeof(x));
    }

    template <class X>
    static void read(std::ifstream &f, X &x)
    {
		f.read((char*)&x, sizeof(x));
    }

    static void save_features(std::ofstream &f, const Features &fs)
    {
        if (!fs.size())
            throw std::runtime_error("saving empty features");

        write(f, fs.size());
        write(f, fs.norm);
        for (auto &x : fs)
        {
            write(f, x.first);
            write(f, x.second);
        }
    }

    static Features load_features(std::ifstream &f)
    {
        Features fs;
        size_t size = 0;
        read(f, size);
        read(f, fs.norm);
        if (!size || !fs.norm)
            throw std::runtime_error("loading empty or 0 norm features");

        for (int k = 0; k < int(size); k++)
        {
            LowerLevelFeatureIndex i = 0;
            read(f, i);

            LowerLevelFeatureCount c = 0;
            read(f, c);
            if (!c)
                throw std::runtime_error("load_features:: lower level feature count is 0");
            fs[i] = c;
        }

        return fs;
    }

public:

    static const int VERSION = 2;

    static void save_tree(const Tree &tree, const std::string &filename, int levels)
    {
        std::ofstream f(filename, std::ios::out | std::ios::binary);

        write(f, VERSION);

        if (!tree.empty() && levels > 0)
            save_node(f, tree, tree.head_id, levels, 0);
    }

    static void save_node(std::ofstream &f, const Tree &tree, Tree::node_id node, int levels, int level)
    {
        auto &n = tree.get_node(node);
        save_features(f, n.features);
        bool final = n.final() || (level >= levels - 1);
        write(f, final);
        if (!final)
        {
            save_node(f, tree, n.a, levels, level + 1);
            save_node(f, tree, n.b, levels, level + 1);
        }
    }

    static void load_tree(Tree &tree, const std::string &filename)
    {
        tree.nodes.clear();
        tree.nodes.reserve(1000000); // todo: think
        tree.head_id = Tree::INVALID_NODE;

        std::ifstream f(filename, std::ios::in | std::ios::binary);
        if (!f.good())
            throw std::runtime_error("cannot load tree");

        {
            int version = 0;
            read(f, version);
            if (version != VERSION)
                throw std::runtime_error("unsupported tree version");
        }

        if (!f.eof())
            tree.head_id = load_node(f, tree, Tree::INVALID_NODE);
    }
    
    static Tree::node_id load_node(std::ifstream &f, Tree &tree, Tree::node_id parent)
    {
        auto features = load_features(f);
        if (features.empty())
            throw std::runtime_error("empty features on load_node");

        bool final = false;
        read(f, final);
        if (final)
        {
            tree.nodes.push_back(Tree::Node(tree.nodes.size(), parent, features));
            return tree.nodes.size() - 1;
        }
        else
        {
            tree.nodes.push_back(Tree::Node(tree.nodes.size(), parent, features, &tree));

            Tree::node_id id = tree.nodes.size() - 1;
            auto a = load_node(f, tree, id);
            auto b = load_node(f, tree, id);

            tree.get_node(id).a = a;
            tree.get_node(id).b = b;

            return id;
        }
    }

    static bool file_exists(const std::string &filename)
    {
        std::ifstream f(filename);
        return f.good();
    }
};

#endif
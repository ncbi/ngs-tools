#ifndef FEATURE_TREE_H_INCLUDED
#define FEATURE_TREE_H_INCLUDED

#include <iostream>
#include <vector>
#include "features.h"

struct Tree
{
    typedef long int node_id;
    static const node_id INVALID_NODE = -1;

    struct Node
    {
        node_id id, parent, a, b;
        Features features;
        const Tree *tree;

        Node(node_id id, node_id parent, node_id a, node_id b, const Tree *tree) : id(id), parent(parent), a(a), b(b), tree(tree) { recalculate_features(); }
        Node(node_id id, node_id parent, const Features &features, const Tree *tree = nullptr) : id(id), parent(parent), a(INVALID_NODE), b(INVALID_NODE), features(features), tree(tree){}

        bool final() const
        {
            return a == INVALID_NODE;
        }

        void recalculate_features()
        {
            if (!tree)
            {
                if (!final())
                    throw std::runtime_error("no tree set for non final node");
                return;
            }

            features = tree->get_node(a).features;
            FeatureOperations::add(features, tree->get_node(b).features);
        }
    };

    typedef std::vector<Node> Nodes;
    Nodes nodes;
    node_id head_id;

    Tree() : head_id(INVALID_NODE)
    {
        static_assert(sizeof(node_id) >= sizeof(size_t), "sizeof(node_id) >= sizeof(size_t)");
        nodes.reserve(10000000);
    }

    bool empty() const
    {
        return nodes.empty();
    }

    bool valid_node(node_id id) const
    {
        return id >=0 && id < nodes.size();
    }

    void set_head(const Features &features)
    {
        if (!nodes.empty())
            throw std::runtime_error("Tree::set_head !nodes.empty");

        head_id = add_leaf_node(features);
    }

    node_id add_leaf_node(const Features &features)
    {
        nodes.push_back(Node(nodes.size(), INVALID_NODE, features));
        return nodes.rbegin()->id;
    }

    const Node &head() const
    {
        return get_node(head_id);
    }

    const Node &get_node(node_id id) const
    {
        if (!valid_node(id))
            throw std::runtime_error("trying to get invalid node");

        return nodes[id];
    }

    int get_height_level(node_id id) const
    {
        int level = -1;
        for_node_and_parents_do(id, [&](const Node &node) { level++; });
        return level;
    }

    Node &get_node(node_id id)
    {
        if (!valid_node(id))
            throw std::runtime_error("trying to get invalid node");

        return nodes[id];
    }

    void merge_nodes(node_id target_id, const Features &features)
    {
        Node &target = get_node(target_id);
        if (!target.final())
            throw std::runtime_error("merge to non final node");

        FeatureOperations::add(target.features, features);
//		for_node_and_parents_do(target.parent, [&](Node &node) { node.recalculate_features(); });
        recalculate_features_for_node_and_parents(target.parent);
    }

    node_id add_group(node_id target, const Features &features, node_id parent)
    {   
        auto leaf_node = add_leaf_node(features);
        nodes.push_back(Node(nodes.size(), parent, target, leaf_node, this));
        auto new_group = nodes.rbegin()->id;
        
        get_node(target).parent = new_group;
        get_node(leaf_node).parent = new_group;
        
        return new_group;
    }

    void insert_neighbour(node_id target_id, const Features &features)
    {
//        std::cout << "adding neighbout to " << target_id << std::endl;
        auto old_target_parent = get_node(target_id).parent;
        auto new_group = add_group(target_id, features, old_target_parent);

        if (target_id == head_id)
        {
            head_id = new_group;
            return;
        }

        replace_child(old_target_parent, target_id, new_group);
        check_balance_for_node_and_parents(new_group);
    }

    void replace_child(node_id parent_id, node_id was, node_id now)
    {
        auto &parent = get_node(parent_id);
        if (parent.a == was)
            parent.a = now;
        else if (parent.b == was)
            parent.b = now;
        else
            throw std::runtime_error("replace_child:: old node not found");
        
		recalculate_features_for_node_and_parents(parent_id);
    }

    void recalculate_features_for_node_and_parents(node_id node)
    {
		for_node_and_parents_do(node, [&](Node &node) { node.recalculate_features(); });
    }

    void check_balance_for_node_and_parents(node_id node)
    {
		for_node_and_parents_do(node, [&](Node &node) { check_balance(node.id); });
    }

    static double weight(Node &n)
    {
        return n.features.norm;
    }

    void check_balance(node_id node)
    {
//        std::cout << "checking balance for " << node << std::endl;
        auto &n = get_node(node);
        if (n.final())
            throw std::runtime_error("checking balance for final node");

        if (n.id == head_id)
            return;

        auto &sibling = get_sibling_node(n);
        auto &a = get_node(n.a);
        auto &b = get_node(n.b);

//        if (node == 4)
  //          std::cout << weight(a) << " " << weight(b) << " " << weight(sibling) << " sibling: " << sibling.id << std::endl;

        if (weight(a) > weight(b))
        {

            if (weight(sibling) < weight(a))
                rebalance(a, sibling);
        }
        else
        {
            if (weight(sibling) < weight(b))
                rebalance(b, sibling);
        }
    }

    Node &get_sibling_node(Node &x)
    {
        auto &p = get_node(x.parent);
        if (p.a == x.id)
            return get_node(p.b);
        else
        {
            if (p.b != x.id)
                throw std::runtime_error("get_sibling_node p.b != x.id");
            return get_node(p.a);
        }
    }

    void rebalance(Node &a, Node &sibling)
    {
//        std::cout << "rebalancing! " << a.id << " " << sibling.id << std::endl;
        auto &b = get_sibling_node(a);
        auto &n = get_node(a.parent);
        auto &p = get_node(n.parent);
        if (get_sibling_node(n).id != sibling.id)
            throw std::runtime_error("rebalance:: get_sibling(n).id != sibling.id");

        p.a = a.id;
        p.b = n.id;

        a.parent = p.id;
        n.parent = p.id;

        n.a = b.id;
        n.b = sibling.id;

        b.parent = n.id;
        sibling.parent = n.id;

        n.recalculate_features();
        p.recalculate_features();
    }

    template <class Lambda>
	void for_node_and_parents_do(node_id node, Lambda &&lambda)
    {
        while (valid_node(node))
        {
            auto &n = get_node(node);
            lambda(n);
            node = n.parent;
        }
    }

    template <class Lambda>
	void for_node_and_parents_do(node_id node, Lambda &&lambda) const
    {
        while (valid_node(node))
        {
            auto &n = get_node(node);
            lambda(n);
            node = n.parent;
        }
    }
};

/*
template <class IDCast = LowerLevelFeatureIndex>
static void print_tree(const Tree &tree, Tree::node_id id, int indent = 0)
{
    const auto &node = tree.get_node(id);
    //if (indent > 5)
    //    return;

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

        print_tree(tree, node.a, indent + 1);
        print_tree(tree, node.b, indent + 1);
    }
}
*/

#endif
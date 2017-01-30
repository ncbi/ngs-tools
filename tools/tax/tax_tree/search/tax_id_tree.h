#ifndef TAX_ID_TREE_H_INCLUDED
#define TAX_ID_TREE_H_INCLUDED

#include <set>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <algorithm>

typedef unsigned int tax_id_t;

struct TaxIdTree
{
	static const tax_id_t ROOT = 1;

	struct Node
	{
		tax_id_t tax_id, parent_tax_id;
		std::set<tax_id_t> subids;

		Node(tax_id_t tax_id, tax_id_t parent_tax_id) : tax_id(tax_id), parent_tax_id(parent_tax_id){}
	};

	typedef std::map<tax_id_t, Node*> Nodes;
//	typedef std::map<tax_id_t, std::unique_ptr<Node> > Nodes;
	Nodes nodes;

	tax_id_t consensus_of(tax_id_t tax_a, tax_id_t tax_b) const
	{
//		std::cerr << "consensus of " << tax_a << " and " << tax_b << std::endl;
		if (tax_a == ROOT || tax_b == ROOT)
			return ROOT;

//		if (tax_a == tax_b)
//			return tax_a;

		if (a_sub_b(tax_a, tax_b))
			return tax_b;

		if (a_sub_b(tax_b, tax_a))
			return tax_a;

		return consensus_of(get_parent_id(tax_a), get_parent_id(tax_b));
	}

	bool a_sub_b(tax_id_t tax_a, tax_id_t tax_b) const
	{
		if (tax_a == tax_b || tax_b == ROOT)
			return true;

		auto node_ptr = get_node_ptr(tax_b);
		return node_ptr->subids.find(tax_a) != node_ptr->subids.end();
	}

	Node *get_node_ptr(tax_id_t tax_id) const
	{
		auto it = nodes.find(tax_id);
		if (it == nodes.end())
		{
			auto message = std::string("no such tax_id as ") + std::to_string(tax_id);
			std::cerr << message << std::endl;
			throw message;
		}

		return it->second;
	}

	tax_id_t get_parent_id(tax_id_t tax_id) const
	{
		return get_node_ptr(tax_id)->parent_tax_id;
	}
};

struct TaxIdTreeLoader
{
	static void load_tax_id_tree(TaxIdTree &tax_id_tree, const std::string &filename)
	{
		std::ifstream f(filename);
		if (f.fail() || f.eof())
			throw std::string("cannot open file ") + filename;

		while (!f.eof())
		{
			tax_id_t x, parent;
			f >> x;
			f >> parent;

			if (!x || !parent)
				throw std::string("bad tax id: ") + std::to_string(x);

			tax_id_tree.nodes[x] = new TaxIdTree::Node(x, parent); // yes, we should free this memory. later...
		}

		calculate_subids(tax_id_tree);
	}

#if 0
	static void	add_subids_recursive(TaxIdTree &tax_id_tree, tax_id_t to, tax_id_t what_to_add, std::set<tax_id_t> &what_to_add_more)
	{
		if (to == TaxIdTree::ROOT)
			return;

		auto to_ptr = tax_id_tree.get_node_ptr(to);
		auto before = to_ptr->subids.size();
		{
			to_ptr->subids.insert(what_to_add);
			to_ptr->subids.insert(what_to_add_more.begin(), what_to_add_more.end());
		}

		if (to_ptr->subids.size() != before)
			add_subids(tax_id_tree, to_ptr->parent_tax_id, to_ptr->tax_id, to_ptr->subids);
	}

	static void calculate_subids_recursive(TaxIdTree &tax_id_tree) // todo: optimize
	{
		for (auto &n : tax_id_tree.nodes)
		{
			auto tax_id = n.first;
			auto *node = n.second;
			add_subids(tax_id_tree, node->parent_tax_id, tax_id, node->subids);
		}
	}
#else

	static void calculate_subids(TaxIdTree &tax_id_tree)
	{
		std::set<TaxIdTree::Node*> nodes;
		for (auto &n : tax_id_tree.nodes)
			nodes.insert(n.second);

		calculate_subids(tax_id_tree, nodes);
	}

	static void calculate_subids(TaxIdTree &tax_id_tree, const std::set<TaxIdTree::Node*> &nodes)
	{
//		std::cout << "calc subids for " << nodes.size() << std::endl;

		std::set<TaxIdTree::Node*> parents;
		for (auto n : nodes)
			if (n->parent_tax_id != TaxIdTree::ROOT)
				parents.insert(tax_id_tree.get_node_ptr(n->parent_tax_id));

		//std::set<TaxIdTree::Node*> diff;
		//auto it = set_difference(nodes.begin(), nodes.end(), parents.begin(), parents.end(), std::inserter(diff, diff.end()));
		//for (; it!= diff.end(); it++)
		for (auto n : nodes)
			if (parents.find(n) == parents.end())
			{
	//			TaxIdTree::Node *n = *it;
				add_subids(tax_id_tree, n->parent_tax_id, n->tax_id, n->subids);
			}

		if (!parents.empty())
			calculate_subids(tax_id_tree, parents);
	}

	static void	add_subids(TaxIdTree &tax_id_tree, tax_id_t to, tax_id_t what_to_add, const std::set<tax_id_t> &what_to_add_more)
	{
		if (to == TaxIdTree::ROOT)
			return;

		auto to_ptr = tax_id_tree.get_node_ptr(to);

		to_ptr->subids.insert(what_to_add);
		to_ptr->subids.insert(what_to_add_more.begin(), what_to_add_more.end());
	}

#endif

};

#endif

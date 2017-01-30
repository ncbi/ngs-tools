"use strict"

var tax_tree = {};

function tree_transform(oData) 
{
//	return oData;
	var tree = tax_tree.to_tree(oData);
	tax_tree.sort_nodes(tree);
	tax_tree.remove_identification_errors(tree);
	tax_tree.remove_percent_duplicates(tree);
	tax_tree.collapse_unclear(tree);
	tax_tree.fix_zero_percent_message(tree);
        return tax_tree.from_tree(tree);
}

(function (tax_tree) {

	tax_tree.to_tree = function (oData)
	{
        var tree = this.tree_node("root", null, {});
        var node_by_name = {}

        for (var i = 0; i < oData.length; i++)
                this.add_to_tree(tree, oData[i], node_by_name);

//        console.info(tree);
        return tree;
	},

	tax_tree.is_root = function (node)
	{
		return node.parent == null && node.name == "root";
	},

	tax_tree.from_tree = function (tree)
	{
		var data = [];
		var self = this;
		this.for_every_tree_node(tree, function (node) {
			if (node != tree)
				data.push(self.node_to_data(tree, node));
			return true;
		});

		return data;
	},

	tax_tree.for_every_tree_node = function (tree, f)
	{
		if (!tree)
			return;

		if (!this.is_root(tree))
			if (!f(tree))
				return;

		for (var i = 0; i < tree.childs.length; i++)
			this.for_every_tree_node(tree.childs[i], f);
	},

	tax_tree.for_every_node_parent = function (node, f)
	{
		while (node)
		{
			node = node.parent;
			if (node)
				f(node);
		}
	},

	tax_tree.tree_find_node_by_name = function (node_by_name, name)  // todo: remove checks?
	{
		if (name in node_by_name) // todo: remove check?
			return node_by_name[name];

		return null;
	},

	tax_tree.node_to_data = function node_to_data(tree, node)
	{
		var x = {n : node.name, p: node.parent != tree ? node.parent.name : "0", d: node.display}
		if (node.collapsed)
			x.collapsed = node.collapsed;
			
		if (node.childs.length == 0)
			x.c = 0;

		return x;
	},

	tax_tree.tree_node = function (tname, tparent, tdisplay)
	{
		return { name : tname, parent : tparent, display : tdisplay, childs : [] };
	}

	tax_tree.add_to_tree = function (tree, node_data, node_by_name)
	{
		if (node_data.n === undefined) // our data bug
			return;
			
		var self = this;
		var node = this.tree_node(node_data.n, self.tree_find_node_by_name(node_by_name, node_data.p), node_data.d);
		if (!node.parent)
			node.parent = tree;

		node.parent.childs.push(node);

		if (node.name in node_by_name)
			throw "tree nodes have duplicate names";

		node_by_name[node.name] = node;
	},

	tax_tree.sort_nodes = function (tree)
	{
		this.for_every_tree_node(tree, function (node) {
			node.childs.sort( function (a, b) { return b.display.percent - a.display.percent; } );
			return true;
		});
		
		tree.childs.sort( function (a, b) { return b.display.percent - a.display.percent; } );
	},

	tax_tree.promote_child = function (tree, node)
	{
		var parent = node.parent;
		var child = node.childs[0]
		if (!parent) // || node.childs[0].length != 1)
			return null;

		parent.childs = node.childs;
		child.parent = parent;	
		return child;
	},

	tax_tree.remove_node = function remove_node(tree, node)
	{
		var parent = node.parent;
		var node_i = parent.childs.indexOf(node)
		if (node_i < 0)
			return; //throw "cannot remove node";

		parent.childs.splice(node_i, 1);
		node.childs = []
	},

	tax_tree.percent_almost_similar = function (node_percent, child_percent)
	{
	 	return Math.abs(node_percent - child_percent) <= node_percent/20.0;
	}

	tax_tree.remove_percent_duplicates = function (tree)
	{
		var self = this;
		this.for_every_tree_node(tree, function (node) {
			while (node.parent && 
				node.parent.childs.length == 1 && 
				node.childs.length == 1 && 
				self.percent_almost_similar(node.display.percent, node.childs[0].display.percent) && //				node.childs[0].display.percent == node.display.percent && 
				(node = self.promote_child(tree, node)));

			return true;
		});
	},

	tax_tree.is_eukaryota = function (from_node)
	{
		if (from_node.display.name == "Eukaryota")
			return true;

		var euk_found = false;
		this.for_every_node_parent(from_node, function (node) {
			if (node.display.name == "Eukaryota")
				euk_found = true;
		});
			
		return euk_found;	
	},

	tax_tree.is_human = function (from_node)
	{
		var human_found = false;
		this.for_every_tree_node(from_node, function (node) {
			if (node.display.name == "Homo sapiens")
				human_found = true;

			return !human_found;
		});

		return human_found;
	},

	tax_tree.get_top_percent = function (childs)
	{
		if (!childs)
			return 0;
		var top = 0;
		for (var i = 0; i < childs.length; i++)
			if (childs[i].display.percent > top)
				top = childs[i].display.percent;

		return top;
	},

	tax_tree.node_level = function (node)
	{
		var level = -1;
		while (node)
		{
			node = node.parent;
			level++;
		}

		return level;
	},

	tax_tree.collapse_unclear = function (tree)
	{
		var MIN_COLLAPSE_PERCENT = 5;
		var MAX_CHILDS = 5;
		var MAX_CHILDS_DEEP = 3;
		var DEEP_NODE_LEVEL = 4;
		var self = this;
		this.for_every_tree_node(tree, function (node) {
			if (node.childs.length > 0 && node.display.percent < MIN_COLLAPSE_PERCENT)
			{
				node.collapsed = true;
				return false;
			}

			if (node.childs.length >= MAX_CHILDS || (node.childs.length >= MAX_CHILDS_DEEP && self.node_level(node) >= DEEP_NODE_LEVEL) )
				node.collapsed = true;

			return true;
		});
	},


	tax_tree.remove_identification_errors = function (tree)
	{
		var MIN_PERCENT = 1.0;
		var MIN_NODE_LEVEL = 5;
		var MIN_TOP_PART = 20.0;

		var self = this;
		var to_remove = []
		this.for_every_tree_node(tree, function (node) {
			if (node.parent && 
				node.parent.childs.length > 1 && 
				( (node.display.percent < MIN_PERCENT) || (node.display.percent <= self.get_top_percent(node.parent.childs)/MIN_TOP_PART)) &&// && node.childs.length == 1 && node.childs[0].display.percent == node.display.percent && (node = promote_child(tree, node)));
				self.node_level(node) > MIN_NODE_LEVEL && 
				self.is_eukaryota(node) && 
				!self.is_human(node))
			{
//				console.info("to remove");
//				console.info(node);
//				self.remove_node(tree, node)
				to_remove.push(node);
			}
			
			return true;
		});
		
		for (var i = 0; i < to_remove.length; i++)
			this.remove_node(tree, to_remove[i]);
	},

	tax_tree.fix_zero_percent_message = function (tree)
	{
		this.for_every_tree_node(tree, function (node) {
			if (node.display.percent == 0)
				node.display.percent = "< " + 0.01;
			return true;
		});
	}
}(tax_tree));


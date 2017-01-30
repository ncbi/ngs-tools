#!/usr/bin/python
import sys
import os

class Leaf:
	def __init__(self, name):
		self.name = name
		self.count = 1

	def final(self):
		return True

class Group:
	def __init__(self, name, childs, distance, count):
		self.name = name
		self.childs = childs
		self.distance = distance
		self.count = count

	def final(self):
		return False


def load_base_nodes(filename):
	nodes = []
	f = open(filename)
	for line in f:
		name = line.rstrip().split('\t')[-1]
		nodes.append(Leaf(name[7:])) #todo: remove hardcode

	return nodes

def common_name(s):
	return ""


def common_name_pair(sa, sb):
	def _iter():
		for a, b in zip(sa, sb):
			if a == b:
				yield a
			else:
				return

	return ''.join(_iter())

def add_merged_nodes(nodes, filename):
	f = open(filename)
	for line in f:
		line = line.split()
		if len(line) < 4:
			break

		childs = line[2 : -1]
		if len(childs) != int(line[1]):
			raise "bad format"
		childs = [int(x) for x in childs]

		distance = float(line[-1])
		name = common_name([nodes[a].name for a in childs])
		count = sum([nodes[a].count for a in childs]) #nodes[a].count + nodes[b].count
		nodes.append(Group(name, childs, distance, count))

	return nodes

def print_node_xml(nodes, number, indent = 0):
	node = nodes[number]
#<taxon id="10" name="name" total_count="10" />
#</taxon>
	for i in xrange(indent):
		print '\t',

	if node.final():
		print '<taxon id="' + str(number) + '" name="' + node.name + '" />'
		return

	print '<taxon id="' + str(number) + '" name="' + node.name + '" total_count="' + str(node.count) + '" distance="' + str(node.distance) + '" >'
	for a in node.childs:
		print_node_xml(nodes, a, indent + 1)

	for i in xrange(indent):
		print '\t',

	print '</taxon>'

def print_node_newick(nodes, number):
	node = nodes[number]

	if node.final():
		print node.name, 
		return

	print '(',
	print_node_newick(nodes, node.a)
	print ',',
	print_node_newick(nodes, node.b)
	print ')',

def print_tree(nodes):
	print_node_xml(nodes, len(nodes) - 1)
#	print_node_newick(nodes, len(nodes) - 1)

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 3:
		print "need <files.list> <merge list file>" 
		return

	nodes = load_base_nodes(sys.argv[1])
	add_merged_nodes(nodes, sys.argv[2])
	print_tree(nodes)
	
main()


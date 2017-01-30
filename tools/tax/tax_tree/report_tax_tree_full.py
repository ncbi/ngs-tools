#!/opt/python-2.7/bin/python
#!/usr/bin/python
import os
import sys

class Root:
	def __init__(self):
		self.nodes = dict()
		self.weight = 0

def build_tree( lineages_count ):
	roots = dict()
	return add_tree_level(roots, lineages_count, 0)

def add_tree_level(roots, lineages_count, level):
	sum_count = 0
	for lineage, count in lineages_count:
		if level >= len(lineage):
			continue
		line = lineage[level]
		if not line in roots:
			roots[line] = Root()

		roots[line].weight += count
		sum_count += count

	for line in roots:
		same_lineage_root = [ (lineage, count) for (lineage, count) in lineages_count if len(lineage) > level and lineage[level] == line ]
		add_tree_level(roots[line].nodes, same_lineage_root, level + 1)

	return roots, sum_count


def one_and_equal_element(roots, count):
	if not roots or len(roots) != 1:
		return False

	roots = [roots[line].weight for line in roots]
	return roots[0] == count

def print_tree_wrong(roots, level = 0):
	roots = [(line, roots[line].nodes, roots[line].weight) for line in roots]
	roots.sort(key = lambda(line, subtree, count) : - count)
	for line, subtree, count in roots:
		if not one_and_equal_element(subtree, count):
			for i in range(level):
				print '\t',
			print line, count
		print_tree(subtree, level + 1)
	
def print_tree_full(roots, level = 0):
	roots = [(line, roots[line].nodes, roots[line].weight) for line in roots]
	roots.sort(key = lambda(line, subtree, count) : - count)
	for line, subtree, count in roots:
		for i in range(level):
			print '\t',
		print line, count
		print_tree(subtree, level + 1)
	
def print_tree(roots, level = 0):
	roots = [(line, roots[line].nodes, roots[line].weight) for line in roots]
	roots.sort(key = lambda(line, subtree, count) : - count)
	if not roots:
		return

	top_count = roots[0][2]

	for line, subtree, count in roots:
		for i in range(level):
			print '\t',

		print line, count
		print_tree(subtree, level + 1)
	
	
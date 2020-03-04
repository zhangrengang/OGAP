import sys
import re
from Bio import Phylo
def convert_ete3_ncbi_tree_to_species_tree(ncbi_tree, species_tree, trefmt='newick'):
	tree = Phylo.read(ncbi_tree, trefmt)
	for clade in tree.get_terminals() + tree.get_nonterminals():
		try:
			name = re.compile(r':name=(.*?) \- \d+:rank').search(clade.comment).groups()[0]
			name = name.replace(' ', '_')
			clade.name = name
			clade.comment = None
			clade.confidence = None
			clade.branch_length = 0.1
		except TypeError: pass
	Phylo.write(tree, species_tree, trefmt)

if __name__ == '__main__':
	convert_ete3_ncbi_tree_to_species_tree(ncbi_tree=sys.argv[1], species_tree=sys.argv[2])

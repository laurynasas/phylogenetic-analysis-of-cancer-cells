from numpy import array, zeros, append, round, int8
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyTreeConstructor, ParsimonyScorer, NNITreeSearcher
from Bio.Phylo import draw_ascii,draw_graphviz,to_networkx
from Bio import AlignIO
from Bio import Phylo
import operator
import networkx as nx
import matplotlib.pyplot as plt

if __name__ == '__main__':
    aln = AlignIO.read('./data/genotype.phy', 'phylip')
    print aln
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    print dm
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    print "Neighbour joining:"
    print(tree)
    draw_ascii(tree)
    # Phylo.draw_graphviz(tree)


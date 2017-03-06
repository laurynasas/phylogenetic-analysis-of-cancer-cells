from numpy import array, zeros, append, round, int8
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyTreeConstructor, ParsimonyScorer, NNITreeSearcher
from Bio.Phylo import draw_ascii,draw_graphviz,to_networkx
from Bio import AlignIO
import operator
import networkx as nx


def get_genotypes_from_clusters(cluster_dict, vector_size):
    genotype = []
    for key, values in cluster_dict.items():
        genotype_sum = array(zeros(vector_size))
        for sample in values:
            genotype_sum += array(sample)
        # print "-->", genotype_sum / float(len(values))
        genotype += [(array(round(genotype_sum / float(len(values))), dtype=int8)).tolist()]
    dict = {}
    print "^^^",len(genotype),genotype
    for el in genotype:
        dict[sum(el)] = el
    sorted_x = sorted(dict.items(), key=operator.itemgetter(0))

    return [el[1] for el in sorted_x]

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


    print "--------------------------"
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor_parsimony = ParsimonyTreeConstructor(searcher)
    pars_tree = constructor_parsimony.build_tree(aln)
    print(pars_tree)
    networkx_tree = to_networkx(pars_tree)
    draw_graphviz(pars_tree)


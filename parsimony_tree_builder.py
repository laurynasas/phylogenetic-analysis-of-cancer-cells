import operator

import numpy as np
from Bio import AlignIO
from Bio.Phylo import to_networkx
from Bio.Phylo.TreeConstruction import ParsimonyTreeConstructor, ParsimonyScorer, NNITreeSearcher
from numpy import array, zeros, round, int8


def get_genotypes_from_clusters(cluster_dict, vector_size):
    genotype = []
    for key, values in cluster_dict.items():
        genotype_sum = array(zeros(vector_size))
        for sample in values:
            genotype_sum += array(sample.split(','), dtype=int8)
        # print "-->", genotype_sum / float(len(values))
        print genotype_sum
        genotype += [(array(round(genotype_sum / float(len(values))), dtype=int8)).tolist()]
    dict = {}
    print "^^^", len(genotype), genotype
    for el in genotype:
        dict[sum(el)] = el
    sorted_x = sorted(dict.items(), key=operator.itemgetter(0))

    return [el[1] for el in sorted_x]


def get_distance_between_genotypes(x, y):
    distance = 0
    for i in xrange(len(x)):
        if x[i] != y[i]:
            distance += 1
    return distance


def prepare_raw_genotypes_for_tree_reconstruction(full_true_genotype, full_predicted_geotype,
                                                  allowed_hamming_distance=None):
    predicted_genotypes = []
    full_predicted_geotype = np.array(full_predicted_geotype)
    for pred_gen in full_predicted_geotype:
        all_distances = []
        for true_gen in full_true_genotype:
            all_distances += [get_distance_between_genotypes(true_gen, pred_gen)]
        all_distances = np.array(all_distances)
        if allowed_hamming_distance != None:
            possible_distance_set_indices = np.where(all_distances <= allowed_hamming_distance)
            if len(possible_distance_set_indices[0]) == 0:
                continue
        predicted_genotypes.append(pred_gen)
    return predicted_genotypes


if __name__ == '__main__':
    aln = AlignIO.read('./data/genotype.phy', 'phylip')

    print "--------------------------"
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor_parsimony = ParsimonyTreeConstructor(searcher)
    pars_tree = constructor_parsimony.build_tree(aln)
    print(pars_tree)
    networkx_tree = to_networkx(pars_tree)
    print networkx_tree.nodes()
    # draw_ascii(pars_tree)
    # Phylo.draw_graphviz(pars_tree)
    # networkx_tree = networkx_tree.to_directed()
    # networkx_tree = nx.dodecahedral_graph()
    # nx.draw_networkx(networkx_tree)
    # nx.draw(networkx_tree, pos=nx.spring_layout(networkx_tree))
    # plt.show()
    # draw_graphviz(pars_tree)

    full_true_genotype = [[0, 0, 0], [0, 0, 1]]
    full_pred_genotype = [[0, 1, 0]]
    print prepare_raw_genotypes_for_tree_reconstruction(full_true_genotype, full_pred_genotype,
                                                        allowed_hamming_distance=0)

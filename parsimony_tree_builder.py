from Bio import AlignIO
from Bio.Phylo import to_networkx, draw_ascii
from Bio.Phylo.TreeConstruction import ParsimonyTreeConstructor, ParsimonyScorer, NNITreeSearcher
from numpy import int8
from Bio import Phylo
from k_medoid_clustering import *
import operator
import pylab as plt
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyTreeConstructor, ParsimonyScorer, NNITreeSearcher
import networkx as nx
from Bio.Alphabet import *



def get_genotypes_from_clusters(cluster_dict, vector_size ,delimiter = ","):
    genotype = []
    for key, values in cluster_dict.items():
        genotype_sum = array(zeros(vector_size))
        for sample in values:
            genotype_sum += array(sample.split(delimiter), dtype=int8)
        # print "-->", genotype_sum / float(len(values))
        # print genotype_sum
        genotype += [(array(np.round(genotype_sum / float(len(values))), dtype=int8)).tolist()]
    dict = {}
    # print "^^^", len(genotype), genotype
    for el in genotype:
        dict[sum(el)] = el
    sorted_x = sorted(dict.items(), key=operator.itemgetter(0))

    return [el[1] for el in sorted_x]


def read_true_genotypes(dir):
    sample_data_file = open(dir, 'rw+')
    lines = sample_data_file.readlines()

    genotypes = []
    for line in lines:
        genotypes.append(map(int, line.split(",")))
    return genotypes


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

def read_original_genotypes(dir):
    target = open(dir)
    data = target.readlines()
    genotypes = []
    for line in data:
        genotypes.append([str(el) for el in line])
    return genotypes

def write_genotype_to_file(dir, ready_gen):
    target = open(dir,  'w+')
    target.write(str(len(ready_gen))+" "+str(len(ready_gen[0])) + "\n")

    for index, el in enumerate(ready_gen):
        target.write(str(index).ljust(10) + ','.join(map(str,el))+"\n")
    target.close()

if __name__ == '__main__':
    true_data_dir = "/home/laurynas/workspace/individual_project/simulated_data/true_genotypes_10_10_0.01_100.txt"
    full_true_genotype = read_true_genotypes(true_data_dir)
    sample_name = "populated_true_genotypes_10_10_0.01_100"
    dir = "/home/laurynas/workspace/individual_project/simulated_data/" + sample_name + ".txt"
    vector_size = 10
    no_clusters = 20
    k_means_instance = k_medoid(dir=dir, number_of_clusters=no_clusters, vector_size=vector_size)
    predicted_clusters = k_means_instance.do_k_means_using_sklearn()
    full_pred_genotype = get_genotypes_from_clusters(predicted_clusters, vector_size, delimiter=", ")
    ready_genotype = prepare_raw_genotypes_for_tree_reconstruction(full_true_genotype, full_pred_genotype,
                                                  allowed_hamming_distance=10)

    sample_name = "predicted_ready_genotypes_10_10_0.01_100"
    dir = "/home/laurynas/workspace/individual_project/simulated_data/" + sample_name + ".phy"
    write_genotype_to_file(dir, ready_genotype)

    sample_name_write = "predicted_SLC_ready_genotypes_10_10_0.01_100"
    sample_name_read = "true_genotypes_10_10_0.01_100"
    read_dir = "/home/laurynas/workspace/individual_project/simulated_data/" + sample_name_read + ".txt"
    write_dir = "/home/laurynas/workspace/individual_project/simulated_data/" + sample_name_write + ".phy"

    true_genotypes = read_true_genotypes(read_dir)
    # write_genotype_to_file(write_dir, [np.array(list_el) for list_el in true_genotypes])

    letters = [0,1]
    alphabet = SingleLetterAlphabet()
    alphabet.letters = letters
    # alphabet = AlphabetEncoder(letters,1)
    aln = AlignIO.read(write_dir, 'phylip', alphabet = alphabet)
    print aln
    # # print aln
    # calculator = DistanceCalculator('identity')
    # dm = calculator.get_distance(aln)
    # print dm
    # constructor = DistanceTreeConstructor(calculator, 'nj')
    # tree = constructor.build_tree(aln)
    # print tree.is_terminal()
    # draw_ascii(tree)
    # aln = AlignIO.read(dir, 'phylip')
    #
    # print "--------------------------"
    dtc = DistanceTreeConstructor(DistanceCalculator("identity"),
                                  "upgma")
    # print dtc.distance_calculator.get_distance()
    nj_tree = dtc.build_tree(aln)
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor_parsimony = ParsimonyTreeConstructor(searcher, starting_tree = nj_tree)
    pars_tree = constructor_parsimony.build_tree(aln)
    # pars_tree.ladderize(False)
    print(pars_tree)
    # print(pars_tree)
    # networkx_tree = to_networkx(pars_tree)
    # # print networkx_tree.nodes()
    # draw_ascii(pars_tree)
    Phylo.draw_graphviz(pars_tree, prog="dot")
    # networkx_tree = networkx_tree.to_directed()
    # networkx_tree = nx.dodecahedral_graph()
    # nx.draw_networkx(networkx_tree)
    # nx.draw(networkx_tree, pos=nx.spring_layout(networkx_tree))
    plt.show()
    # draw_graphviz(pars_tree)

    full_true_genotype = [[0, 0, 0], [0, 0, 1]]
    full_pred_genotype = [[0, 1, 0]]
    # print prepare_raw_genotypes_for_tree_reconstruction(full_true_genotype, full_pred_genotype,
    #                                                     allowed_hamming_distance=0)

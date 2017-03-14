import processing_sample_data as pr
import bernoulli_mixture_model as bm
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

def count_differences(str1, str2):
    count = 0
    str1 = str1.split(",")
    str2 = str2.split(",")
    # print str1
    # print str2
    for el1, el2 in zip(str1, str2):
        if el1 =="NA" or el2 =="NA":
            continue
        if el1 != el2:
            count += 1
    return count


def calc_distance_matrix_from_file(dir):
    sample_data_file = open(dir, 'rw+')
    all_lines = []
    for line in sample_data_file:
        # if len(line)>2:
        all_lines += [line]
    sample_data_file.close()

    matrix = []
    for i in all_lines:
        distance_row = []
        for j in all_lines:
            distance_row += [count_differences(i, j)]
        matrix += [distance_row]
    return matrix


def plot_2D_similarity_matrix(matrix, method="", dataset="", no_cl="", data_type=""):
    from matplotlib import mpl, pyplot

    # make values from -5 to 5, for this example
    # zvals = np.random.rand(100, 100) * 10 - 5

    # make a color map of fixed colors
    cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                         ['blue', 'black', 'red'],
                                                         256)
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    img2 = pyplot.imshow(matrix, interpolation='nearest',
                         cmap=cmap2,
                         origin='lower')

    # make a color bar

    pyplot.colorbar(img2, cmap=cmap2)
    fig = pyplot.figure(1)

    # fig.savefig(
    #     "/home/laurynas/workspace/individual_project/data/clutsering_analysis_images/" + str(data_type) + "/" + str(
    #         method) + "_" + str(dataset) + "_" + str(no_cl) + ".png", bbox_inches='tight')
    # fig.savefig("image2.png")
    pyplot.show()


def get_single_cell_genotypes(clustered_data_dict, delimiter=","):
    genotypes = []
    for key, values in clustered_data_dict.items():
        genotype = [{'0': 0, '1': 0, 'NA': 0} for _ in range(len(values[0].split(delimiter)))]
        print "getting genotypes: ", key, len(values)
        for sample in values:
            for index,el in enumerate(sample.split(delimiter)):
                genotype[index][el] +=1
        genotype_list =[]
        for index, el_dict in enumerate(genotype):
            genotype_list.append(max(el_dict, key=el_dict.get))
        genotypes.append(genotype_list)

    sums = []
    string_genotype = []
    for el in genotypes:
        string_genotype += [",".join(str(e) for e in el)]

    # sorted_sums, string_genotype = (list(t) for t in zip(*sorted(zip(sums, string_genotype))))

    return string_genotype


if __name__ == "__main__":
    directory = "/home/laurynas/workspace/individual_project/data/hou/snv.csv"
    unique_rows, full_or_data, full_info = pr.process_single_cell_data_file(directory)
    bmm = bm.BMM(k_clusters=5, unique_rows=unique_rows, full_data_dict=full_or_data, full_info=full_info,
              number_of_iterations=10)
    data = bmm.do_clustering()
    print data.keys()
    target_dir = "/home/laurynas/workspace/individual_project/data/hou/clusters_5_snv.txt"
    target = open("/home/laurynas/workspace/individual_project/data/hou/clusters_5_snv.txt", 'w+')
    for key in data.keys():
        print key, len(data[key])
        for el in data[key]:
            target.write(el+"\n")

    distance_matrix = calc_distance_matrix_from_file(target_dir)
    # print distance_matrix
    plot_2D_similarity_matrix(distance_matrix)
    plt.show()

    # ready_gen = get_single_cell_genotypes(data)
    # write_dir = "/home/laurynas/workspace/individual_project/data/hou/clusters_5_genotypes_snv.phy"
    # target = open(write_dir, 'w+')
    #
    # target.write(str(len(ready_gen)) + " " + str(len(ready_gen[0])) + "\n")
    #
    # for index, el in enumerate(ready_gen):
    #     target.write(str(index).ljust(10) + el.replace("NA","N") + "\n")
    # target.close()
    #
    #
    #
    # aln = AlignIO.read(write_dir, 'phylip')
    # dtc = DistanceTreeConstructor(DistanceCalculator("identity"),
    #                               "nj")
    # # print dtc.distance_calculator.get_distance()
    # nj_tree = dtc.build_tree(aln)
    # print nj_tree
    # scorer = ParsimonyScorer()
    # searcher = NNITreeSearcher(scorer)
    # constructor_parsimony = ParsimonyTreeConstructor(searcher, starting_tree = nj_tree)
    # pars_tree = constructor_parsimony.build_tree(aln)
    # Phylo.draw_graphviz(pars_tree, prog="dot")
    #
    #
    # plt.show()

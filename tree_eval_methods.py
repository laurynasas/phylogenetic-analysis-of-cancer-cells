from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyTreeConstructor, \
    ParsimonyScorer, NNITreeSearcher
from numpy import int8
from scipy.cluster.hierarchy import fcluster, linkage

from bernoulli_mixture_model import BMM
from binary_tree import BTree
from n_tree import NTree
from k_medoid_clustering import *
from Bio import Phylo
import pylab as plt


class Pipeline:
    def __init__(self, raw_data_dir, true_genotypes_dir, clustering_method, tree_method, number_of_clusters,
                 vector_size, write_adjusted_true_gen_dir, write_predicted_gen_dir, bbm_number_of_iterations=10,
                 starting_tree="upgma"):
        self.clustering_method = clustering_method
        self.raw_data_dir = raw_data_dir
        self._starting_tree = starting_tree
        self.true_genotypes_dir = true_genotypes_dir
        self.no_clusters = number_of_clusters
        self.vector_size = vector_size
        self.tree_method = tree_method
        self.write_adjusted_true_gen_dir = write_adjusted_true_gen_dir
        self.write_predicted_gen_dir = write_predicted_gen_dir
        self.number_of_iterations = bbm_number_of_iterations

    def run_pipe(self):

        self._read_raw_simulated_data_file()

        if self.clustering_method == "slc":
            self._do_slc()
        elif self.clustering_method == "k_means":
            self._do_k_means()
        elif self.clustering_method == "bmm":
            self._do_bmm()

        self._get_non_unique_genotypes()
        self._load_true_genotypes()

        self._equate_dimensions()

        self._save_genotypes()
        self._save_genotypes(true_genotypes=True)

        if self.tree_method == "nj":
            self._build_distance_tree(tree_name="nj")
            self._build_distance_tree(tree_name="nj", true_tree=True)
        elif self.tree_method == "pars":
            self._do_parsimony()
            self._do_parsimony(true_tree=True)

        elif self.tree_method == "upgma":
            self._build_distance_tree(tree_name="upgma")
            self._build_distance_tree(tree_name="upgma", true_tree=True)

        self._get_tree_distances(true_tree=True)
        self._get_tree_distances()

        print np.matrix(self._true_tree_distance_matrix)
        print np.matrix(self._tree_distance_matrix)
        self._get_diff_vector()

        distance = self._calculate_euclidiean_distance()

        return distance

    def _calculate_euclidiean_distance(self):
        sum = 0
        for x, y in zip(self.vector_predicted_diff, self.vector_true_diff):
            sum += abs(x - y) ** 2
        return sum ** (1.0 / 2)

    def _get_diff_vector(self):
        self.vector_predicted_diff = self._get_vector(self._tree_distance_matrix)
        self.vector_true_diff = self._get_vector(self._true_tree_distance_matrix)

    def _get_vector(self, matrix):
        vector = []
        for row in matrix:
            for el in row:
                if el != 0:
                    vector.append(el)
                else:
                    break
        return vector

    def _get_tree_distances(self, true_tree=False):
        if true_tree:
            text_repres = self._true_tree.__str__()
            # if self.tree_method == "pars":
            #     Phylo.draw_graphviz(self._true_tree, prog="dot")
            # else:
            #     Phylo.draw_graphviz(self._true_tree, prog="twopi")
            # plt.show()
        else:
            text_repres = self._tree.__str__()
            # if self.tree_method == "pars":
            #     Phylo.draw_graphviz(self._tree, prog="dot")
            # else:
            #     Phylo.draw_graphviz(self._tree, prog="twopi")
            # plt.show()

        if self.tree_method == "nj" or self.tree_method == "upgma":
            tree = NTree()
        elif self.tree_method == "pars":
            tree = BTree()



        tree.build_tree(lines=text_repres)

        all_nodes = tree.get_all_nodes()
        tree_distance_matrix = []

        for outter_n in all_nodes.values():
            if "Inner" in outter_n.name:
                continue
            row = []
            for inner_n in all_nodes.values():
                if "Inner" in inner_n.name:
                    continue
                if outter_n == inner_n:
                    row.append(0)
                else:
                    row.append(tree.get_distance_between_nodes(outter_n, inner_n))
            tree_distance_matrix.append(row)

        if true_tree:
            self._tree_distance_matrix = tree_distance_matrix
        else:
            self._true_tree_distance_matrix = tree_distance_matrix

    def _do_parsimony(self, true_tree=False):
        if true_tree:
            aln = AlignIO.read(self.write_adjusted_true_gen_dir, 'phylip')
            dtc = DistanceTreeConstructor(DistanceCalculator("identity"),
                                          self._starting_tree)
            starting_tree = dtc.build_tree(aln)
            scorer = ParsimonyScorer()
            searcher = NNITreeSearcher(scorer)
            constructor_parsimony = ParsimonyTreeConstructor(searcher, starting_tree=starting_tree)
            pars_tree = constructor_parsimony.build_tree(aln)
            self._true_tree = pars_tree
        else:
            aln = AlignIO.read(self.write_predicted_gen_dir, 'phylip')
            dtc = DistanceTreeConstructor(DistanceCalculator("identity"),
                                          self._starting_tree)
            starting_tree = dtc.build_tree(aln)
            scorer = ParsimonyScorer()
            searcher = NNITreeSearcher(scorer)
            constructor_parsimony = ParsimonyTreeConstructor(searcher, starting_tree=starting_tree)
            pars_tree = constructor_parsimony.build_tree(aln)
            self._tree = pars_tree

    def _build_distance_tree(self, tree_name, true_tree=False):
        if true_tree:
            aln = AlignIO.read(self.write_adjusted_true_gen_dir, 'phylip')
            calculator = DistanceCalculator('identity')
            constructor = DistanceTreeConstructor(calculator, tree_name)
            tree = constructor.build_tree(aln)
            self._true_tree = tree
        else:
            aln = AlignIO.read(self.write_predicted_gen_dir, 'phylip')
            calculator = DistanceCalculator('identity')
            constructor = DistanceTreeConstructor(calculator, tree_name)
            tree = constructor.build_tree(aln)
            self._tree = tree

    def _read_raw_simulated_data_file(self):
        class DataNode:
            def __init__(self, vector, cluster_label):
                self.vector = vector
                self.cluster_label = cluster_label

        simulated_data_file = open(self.raw_data_dir, 'rw+')
        self.unique_rows = {}
        self.full_data_dict = {}
        self.full_info = []
        for line in simulated_data_file:

            line = line.split(" | ")
            cluster_label = int(line[1])
            vector = line[0][1:-1]
            self.full_info.append(DataNode(vector, cluster_label))

            if self.full_data_dict.get(cluster_label):
                self.full_data_dict[cluster_label] += [map(int, vector.split(","))]
            else:
                self.full_data_dict[cluster_label] = [map(int, vector.split(","))]
            if vector in self.unique_rows.keys():
                self.unique_rows[vector] += 1
            else:
                self.unique_rows[vector] = 1

    def _save_genotypes(self, true_genotypes=False):
        if true_genotypes:
            target = open(self.write_adjusted_true_gen_dir, 'w+')
            ready_gen = self.predicted_genotypes
        else:
            target = open(self.write_predicted_gen_dir, 'w+')
            ready_gen = self._adapted_true_genotypes

        target.write(str(len(ready_gen)) + " " + str(len(ready_gen[0])) + "\n")

        for index, el in enumerate(ready_gen):
            target.write(str(index).ljust(10) + el + "\n")
        target.close()

    def _find_distance_matrix(self):
        number_rows = len(self.unique_rows.keys())
        unique_strings = []
        # unique_strings = unique_rows.keys()

        for key in self.unique_rows.keys():
            for _ in xrange(self.unique_rows[key]):
                unique_strings += [key]

        number_rows = len(unique_strings)

        self.distance_matrix = []
        for i in range(number_rows):
            new_row = []
            for j in range(number_rows):
                new_row += [self._find_number_of_diff_chars_in_string(unique_strings[i], unique_strings[j])]
            self.distance_matrix += [new_row]

        self.distance_matrix = np.matrix(self.distance_matrix)

    def _equate_dimensions(self):
        if len(self.predicted_genotypes) <= len(self.true_genotypes):
            predicted_genotypes_copy = self.predicted_genotypes[:]
            true_genotypes_copy = self.true_genotypes[:]
            self._adapted_true_genotypes = []
            for predicted_gen in predicted_genotypes_copy:
                distances = []
                for true_gen in true_genotypes_copy:
                    distances.append(self._find_number_of_diff_chars_in_string(predicted_gen, true_gen))
                min_index = distances.index(min(distances))
                min_distance_el = true_genotypes_copy[min_index]
                self._adapted_true_genotypes += [min_distance_el]
                true_genotypes_copy.pop(min_index)
        else:
            print "There are more predicted clusters than true clusters!!"

    def _load_true_genotypes(self):
        sample_data_file = open(self.true_genotypes_dir, 'rw+')
        lines = sample_data_file.readlines()

        genotypes = []
        for line in lines:
            genotypes.append(line[:-1])
        self.true_genotypes = genotypes

    def _find_number_of_diff_chars_in_string(self, a, b):
        distance = 0
        for i, j in zip(a, b):
            if i != j:
                distance += 1
        return distance

    def _get_non_unique_genotypes(self, delimiter=","):
        genotype = []
        for key, values in self._clustered_data_dict.items():
            genotype_sum = np.array(np.zeros(self.vector_size))
            for sample in values:
                genotype_sum += np.array(sample.split(delimiter), dtype=np.int8)

            genotype += [(np.array(np.round(genotype_sum / float(len(values))), dtype=np.int8)).tolist()]
        sums = []
        string_genotype = []
        for el in genotype:
            sums += [sum(el)]
            string_genotype += [",".join(str(e) for e in el)]

        sorted_sums, string_genotype = (list(t) for t in zip(*sorted(zip(sums, string_genotype))))

        self.predicted_genotypes = string_genotype

    def _do_slc(self):
        np.matrix(self._find_distance_matrix())
        labels = fcluster(linkage(self.distance_matrix, method='complete'), t=self.no_clusters, criterion='maxclust')
        new_data = {}

        for cluster_label in labels - 1:
            if new_data.get(cluster_label):
                new_data[cluster_label] += [map(int, self.unique_rows.keys()[cluster_label].split(","))]
            else:
                new_data[cluster_label] = [map(int, self.unique_rows.keys()[cluster_label].split(","))]

        new_data_formatetd = {}
        for key, value in new_data.items():
            new_value = [','.join(map(str, el)) for el in value]
            new_data_formatetd[key] = new_value

        self._clustered_data_dict = new_data_formatetd

    def _do_k_means(self):
        predefined_kwargs = {"number_of_clusters": self.no_clusters, "unique_rows": self.unique_rows,
                             "full_data_dict": self.full_data_dict, "full_info": self.full_info,
                             "vector_size": self.vector_size}
        k_means_instance = k_medoid(**predefined_kwargs)

        self._clustered_data_dict = k_means_instance.do_k_means_using_sklearn()

    def _do_bmm(self):
        bmm = BMM(self.no_clusters, self.unique_rows, self.full_data_dict, self.full_info, self.number_of_iterations)
        self._clustered_data_dict = bmm.do_clustering()


raw_data_dir = "/home/laurynas/workspace/individual_project/simulated_data/populated_true_genotypes_10_10_0.01_100.txt"
true_genotype_dir = "/home/laurynas/workspace/individual_project/simulated_data/true_genotypes_10_10_0.01_100.txt"
clustering_method = "bmm"
tree_method = "pars"
number_of_clusters = 10
vector_size = 10
write_adjusted_true_gen_dir = "/home/laurynas/workspace/individual_project/simulated_data/slc_nj_pipe_adjusted_true_gen.phy"
write_predicted_gen_dir = "/home/laurynas/workspace/individual_project/simulated_data/slc_nj_pipe_predicted_gen.phy"

slc_nj_pipe = Pipeline(raw_data_dir, true_genotype_dir, clustering_method, tree_method, number_of_clusters, vector_size,
                       write_adjusted_true_gen_dir, write_predicted_gen_dir)
print slc_nj_pipe.run_pipe()

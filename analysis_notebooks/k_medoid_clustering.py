import numpy as np
from sklearn.cluster import k_means_

import helpers

'''
    The wrapper class for doing k-medoids/k-means clustering.
'''


class KMedoid:
    def __init__(self, *args, **kwargs):
        self.k = kwargs.get('number_of_clusters')
        self.unique_rows = kwargs.get('unique_rows')
        self.full_data_dict = kwargs.get('full_data_dict')
        self.full_info = kwargs.get('full_info')
        self.vector_size = kwargs.get('vector_size')
        if kwargs.get('dir'):
            self.raw_data_dir = kwargs.get('dir')
            self._init_from_dir(kwargs.get('dir'))
        self.sk_learn_instance = None
        self.helper = helpers.Helper()

    '''
        Given the raw data dictionary of simulated deep-sequencing data the methods processes the contents.
    '''

    def _init_from_dir(self, dir):

        unique_rows, full_data_dict, full_info = self.helper.read_simulated_data_file(dir)

        self.unique_rows = unique_rows
        self.full_data_dict = full_data_dict
        self.full_info = full_info

    '''
        Returns predicted labels by my k-medoids implementation.
    '''

    def get_my_predicted_labels(self, C):
        predicted_labels = []
        for key, value in C.items():
            for _ in value:
                predicted_labels.append(key)
        return predicted_labels

    '''
        My implementation of k-medoids
    '''

    def do_k_medoid_my_clustering(self, iterations=100):
        # determine dimensions of distance matrix distance_matrix
        distance_matrix = self.get_distance_matrix()
        number_of_clusters = self.k
        m, n = distance_matrix.shape

        medoids = np.sort(np.random.choice(n, number_of_clusters))
        medoids_copy = np.copy(medoids)
        clusters = {}

        for _ in xrange(iterations):
            cluster_members = np.argmin(distance_matrix[:, medoids], axis=1)
            for cluster_label in range(number_of_clusters):
                clusters[cluster_label] = np.where(cluster_members == cluster_label)[0]
            # update cluster centroids
            for cluster_label in range(number_of_clusters):
                cluster_members = np.mean(distance_matrix[np.ix_(clusters[cluster_label], clusters[cluster_label])],
                                          axis=1)
                if len(cluster_members) != 0:
                    j = np.argmin(cluster_members)
                    medoids_copy[cluster_label] = clusters[cluster_label][j]
            np.sort(medoids_copy)

            # check if it converged
            if np.array_equal(medoids, medoids_copy):
                break

            medoids = np.copy(medoids_copy)
        else:
            # final update of cluster memberships
            cluster_members = np.argmin(distance_matrix[:, medoids], axis=1)
            for cluster_label in range(number_of_clusters):
                clusters[cluster_label] = np.where(cluster_members == cluster_label)[0]

        return medoids, clusters

    '''
        Parses distance matrix to numpy matrix object.
    '''

    def get_distance_matrix(self):
        return np.matrix(self.helper.find_distance_matrix(self.unique_rows))

    '''
        Performs k-means clustering using sklearn library.
    '''

    def do_k_means_using_sklearn(self):
        distance_matrix = self.get_distance_matrix()
        k = k_means_.KMeans(n_clusters=self.k).fit(distance_matrix)
        self.sk_learn_instance = k
        full_dict = {}

        for label, data_object in zip(k.labels_, self.full_info):
            if label not in full_dict.keys():
                full_dict[label] = [data_object.vector]
            else:
                full_dict[label] += [data_object.vector]

        return full_dict

    '''
        Returns cluster labels predicted by sklearn.
    '''

    def get_sklearn_predicted_labels(self):
        return self.sk_learn_instance.labels_

    '''
        Compares two genotypes entry by entry.
    '''

    def compare_genotypes(self, true_genotype, predicted_genotype):
        count = 0
        for gen in predicted_genotype:
            if gen in true_genotype:
                count += 1

        return count / float(len(true_genotype))

    '''
        Loads true genotypes from directory
    '''

    def get_true_genotypes(self, dir):
        simulated_data_file = open(dir, 'rw+')
        data = simulated_data_file.readlines()
        true_genotypes = []
        for line in data:
            true_genotypes.append(map(int, line.split(',')))

        return true_genotypes

    '''
        Finds true labels from the unique rows
    '''

    def get_true_labels(self):

        true_labels = []
        for key in self.unique_rows.keys():
            label = self.helper.get_label_of_cluster(vector=key, full_dict=self.full_data_dict)
            for _ in range(self.unique_rows[key]):
                true_labels.append(label)
        self.true_labels = true_labels

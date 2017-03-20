import numpy as np

import helpers

'''
    Bernoulli Mixture Model clustering class. Please see dissertation class for detailed mathematical description.
'''


class BMM:
    def __init__(self, k_clusters, unique_rows, full_data_dict, full_info, number_of_iterations, distance_matrix=None):
        self.RESTRICTED_CHARS = ["NA\n", "NA"]
        self.unique_rows = unique_rows
        self.k_clusters = k_clusters
        self.full_data_dict = full_data_dict
        self.full_info = full_info
        self.distance_matrix = distance_matrix
        self.number_of_iterations = number_of_iterations
        self.helper = helpers.Helper()
        if not distance_matrix:
            self.distance_matrix = np.matrix(self.helper.find_distance_matrix(unique_rows))

        x_vectors = []
        for key in self.unique_rows.keys():
            rep = self.unique_rows[key]
            x_vector = [str(el) for el in key.split(",")]
            for i in xrange(rep):
                x_vectors.append(x_vector)

        self.x_vectors = np.array(x_vectors)
        self.number_of_vectors = sum(self.unique_rows[key] for key in self.unique_rows.keys())
        self.d_genes = len(self.unique_rows.keys()[0].split(","))

    '''
        Counts the occurrences of elements in a list.
    '''

    def count_occurences(self, list):
        occ = {}
        for el in list:
            try:
                occ[el] += 1
            except:
                occ[el] = 1

        total = sum(occ[key] for key in occ.keys())

        for key in occ.keys():
            occ[key] = occ[key] * 1.0 / total
        return occ

    '''
        Bernoulli probability distribution.
    '''

    def bernoulli_distro(self, p, x):
        return (p ** x) * ((1 - p) ** (1 - x))

    '''
        Total probability that vector x is in cluster k.
    '''

    def prob_x_generated_by_cluster_k(self, x_vector, miu_pixel_distributions, k_cluster):
        product = 1
        for i in xrange(len(x_vector)):
            if x_vector[i] not in self.RESTRICTED_CHARS:
                pwr = int(x_vector[i])
                product *= self.bernoulli_distro(miu_pixel_distributions[k_cluster][i], pwr)
        return product

    '''
        The likelihood function for single mutation vector.
    '''

    def likelihood_of_single_mutation_vector(self, x_vector, pi_proportion_of_mutation_vectors,
                                             miu_pixel_distributions):
        self.k_clusters = len(pi_proportion_of_mutation_vectors)
        likelihood = 0
        for k_cluster in xrange(self.k_clusters):
            likelihood += pi_proportion_of_mutation_vectors[k_cluster] * self.prob_x_generated_by_cluster_k(x_vector,
                                                                                                            miu_pixel_distributions,
                                                                                                            k_cluster)

        return likelihood

    '''
        The overall likelihood function for all data.
    '''

    def likelihood_complete_data(self, pi_proportion_of_vectors, miu_gene_distributions, x_vectors, z_matrix):
        number_of_vectors, d_pixels = x_vectors.shape

        self.k_clusters = len(pi_proportion_of_vectors)

        likelihood = 1
        for n in xrange(number_of_vectors):
            for k in xrange(self.k_clusters):
                likelihood *= np.power(
                    pi_proportion_of_vectors[k] * self.prob_x_generated_by_cluster_k(x_vectors[n],
                                                                                     miu_gene_distributions, k),
                    z_matrix[n, k])
        return likelihood

    '''
        The log of likelihood.
    '''

    def log_likelihood_complete_data(self, likelihood):
        return np.log(likelihood)

    '''
        Used in iterative expectation-maximisation algorithm.
    '''

    def expectation_step(self, z_matrix, pi_proportion_of_vectors, miu_gene_distributions):
        for n in xrange(self.number_of_vectors):
            x_vector = self.x_vectors[n]
            for k in xrange(self.k_clusters):
                nom = pi_proportion_of_vectors[k] * self.prob_x_generated_by_cluster_k(x_vector, miu_gene_distributions,
                                                                                       k)
                den = 0
                for m in xrange(self.k_clusters):
                    den += pi_proportion_of_vectors[m] * self.prob_x_generated_by_cluster_k(x_vector,
                                                                                            miu_gene_distributions, m)

                nom = np.nan_to_num(nom)
                den = np.nan_to_num(den)
                if den == 0:
                    den = 1.0
                z_matrix[n, k] = nom * 1.0 / den

        return z_matrix

    '''
        The result achieved after taking the Lagrangian. See report for explanation.
    '''

    def N_m_effective_number_of_vectors(self, z_matrix, k_cluster):
        total = 0
        for n in xrange(self.number_of_vectors):
            total += z_matrix[n, k_cluster]
        return total

    '''
        Maximization of parameter theta.
    '''

    def maximize_miu(self, miu_gene_distributions, z_matrix):
        for m in xrange(self.k_clusters):
            for j in xrange(self.d_genes):
                N_m = self.N_m_effective_number_of_vectors(z_matrix, m)

                total = 0
                for n in xrange(self.number_of_vectors):
                    if self.x_vectors[n, j] not in self.RESTRICTED_CHARS:
                        total += int(self.x_vectors[n, j]) * z_matrix[n, m]

                miu_gene_distributions[m, j] = total * 1.0 / N_m

        return miu_gene_distributions

    '''
        Maximization of mixing coefficient.
    '''

    def maximize_pi(self, z_matrix, pi_proportion_of_vectors):
        for m in xrange(self.k_clusters):
            N_m = self.N_m_effective_number_of_vectors(z_matrix, m)
            pi_proportion_of_vectors[m] = N_m * 1.0 / self.number_of_vectors

        return pi_proportion_of_vectors

    '''
        Find labels of clustered dictionary.
    '''

    def find_labels(self, z_matrix):
        labels = []

        for i in xrange(z_matrix.shape[0]):
            labels.append(np.argmax(z_matrix[i]))
        return labels

    '''
        Used in gap filling to find possible clusters to take members from.
    '''

    def get_keys_with_more_than_one_el(self, dic):
        keys = []
        for key in dic.keys():
            if len(dic[key]) > 1:
                keys.append(key)
        return keys

    '''
        This is a legacy method. Managed to delete its occurrences, but will leave it in case something crashes.
        It was used to circumvent empty cluster's problem, which was an issue for Silhouette score calculation.
        Since we use sklearn, no problems arose.
    '''

    def fill_the_gaps(self, final_dict, k):

        for i in xrange(k):
            if i not in final_dict.keys():
                if len(final_dict.keys()) == 1:
                    key = final_dict.keys()[0]
                else:
                    possible_keys = self.get_keys_with_more_than_one_el(final_dict)
                    key = np.random.choice(possible_keys)

                final_dict[i] = [final_dict[key].pop(0)]

        return final_dict

    '''
        Main method to combine all of the methods and perform clustering.
    '''

    def do_clustering(self):

        miu_gene_distributions = np.ones((self.k_clusters, self.d_genes))

        for k in xrange(self.k_clusters):
            normalize = 0

            # Init prior distribution
            for i in xrange(self.d_genes):
                miu_gene_distributions[k, i] *= (np.random.random_sample() * 0.5) + 0.25
                normalize += miu_gene_distributions[k, i]

            # Normalize distribution
            for i in xrange(self.d_genes):
                miu_gene_distributions[k, i] /= normalize

        z_matrix = np.zeros((self.number_of_vectors, self.k_clusters), dtype=np.float64)

        z_matrix[:, 0] = np.ones(self.number_of_vectors)
        for i in range(self.number_of_vectors):
            np.random.shuffle(z_matrix[i])

        init_labels = self.find_labels(z_matrix)
        init_dic = {}

        for i in xrange(self.number_of_vectors):
            try:
                init_dic[init_labels[i]].append(self.x_vectors[i])
            except:
                init_dic[init_labels[i]] = []
                init_dic[init_labels[i]].append(self.x_vectors[i])

        pi_proportion_of_vectors = np.sum(z_matrix, axis=0) * 1.0 / self.number_of_vectors

        # Actual iteration of BMM starts
        for i in xrange(self.number_of_iterations):
            # Expectation step
            z_matrix = self.expectation_step(z_matrix, pi_proportion_of_vectors, miu_gene_distributions)

            # Maximization step
            miu_gene_distributions = self.maximize_miu(miu_gene_distributions, z_matrix)
            pi_proportion_of_vectors = self.maximize_pi(z_matrix,
                                                        pi_proportion_of_vectors)

        # Get predicted labels
        final_dict = {}
        labels = self.find_labels(z_matrix)

        for i in xrange(self.number_of_vectors):
            try:
                final_dict[labels[i]].append(",".join(str(el) for el in self.x_vectors[i]))
            except:
                final_dict[labels[i]] = [",".join(str(el) for el in self.x_vectors[i])]


        return final_dict, labels

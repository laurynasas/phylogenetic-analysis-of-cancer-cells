import itertools

from numpy.random import *

from k_medoid_clustering import *
from rand_index import *
import processing_sample_data as pr

class BMM:

    def __init__(self,k_clusters, unique_rows, full_data_dict, full_info, number_of_iterations, distance_matrix=None):
        self.RESTRICTED_CHARS = ["NA\n", "NA"]
        self.unique_rows = unique_rows
        self.k_clusters = k_clusters
        self.full_data_dict = full_data_dict
        self.full_info = full_info
        self.distance_matrix = distance_matrix
        self.number_of_iterations = number_of_iterations

        if not distance_matrix:
            self.distance_matrix = np.matrix(pr.find_distance_matrix(unique_rows))

    def count_occurences(self,list):
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


    def bernoulli_distro(self, p, x):
        return (p ** x) * ((1 - p) ** (1 - x))


    def prob_x_generated_by_cluster_k(self,x_vector, miu_pixel_distributions, k_cluster):
        product = 1
        for i in xrange(len(x_vector)):
            if x_vector[i] not in self.RESTRICTED_CHARS:
                pwr = int(x_vector[i])
                product *= self.bernoulli_distro(miu_pixel_distributions[k_cluster][i], pwr)
        return product


    def likelihood_of_single_training_image(self,x_vector, pi_proportion_of_images, miu_pixel_distributions):
        self.k_clusters = len(pi_proportion_of_images)

        likelihood = 0

        for k_cluster in xrange(self.k_clusters):
            likelihood += pi_proportion_of_images[k_cluster] * self.prob_x_generated_by_cluster_k(x_vector,
                                                                                             miu_pixel_distributions,
                                                                                             k_cluster)

        return likelihood


    def likelihood_complete_data(self,pi_proportion_of_images, miu_pixel_distributions, x_vectors, z_matrix):
        n_images, d_pixels = x_vectors.shape

        self.k_clusters = len(pi_proportion_of_images)

        likelihood = 1
        for n in xrange(n_images):
            for k in xrange(self.k_clusters):
                likelihood *= power(
                    pi_proportion_of_images[k] * self.prob_x_generated_by_cluster_k(x_vectors[n], miu_pixel_distributions, k),
                    z_matrix[n, k])
        return likelihood


    def log_likelihood_complete_data(self,likelihood):
        return log(likelihood)


    def expectation_step(self,z_matrix, pi_proportion_of_images, miu_pixel_distributions):
        for n in xrange(self.n_images):
            x_vector = self.x_vectors[n]
            for k in xrange(self.k_clusters):
                nom = pi_proportion_of_images[k] * self.prob_x_generated_by_cluster_k(x_vector, miu_pixel_distributions, k)
                den = 0
                for m in xrange(self.k_clusters):
                    den += pi_proportion_of_images[m] * self.prob_x_generated_by_cluster_k(x_vector, miu_pixel_distributions, m)

                nom = nan_to_num(nom)
                den = nan_to_num(den)
                if den ==0:
                    den = 1.0
                z_matrix[n, k] = nom * 1.0 / den
        return z_matrix


    def N_m_effective_number_of_images(self,z_matrix, k_cluster):
        total = 0
        for n in xrange(self.n_images):
            total += z_matrix[n, k_cluster]
        return total


    def maximize_miu(self,miu_pixel_distributions, z_matrix, x_vectors):
        for m in xrange(self.k_clusters):
            for j in xrange(self.d_pixels):
                N_m = self.N_m_effective_number_of_images(z_matrix, m)

                total = 0
                for n in xrange(self.n_images):
                    if self.x_vectors[n, j] not in self.RESTRICTED_CHARS:
                        total += int(self.x_vectors[n, j]) * z_matrix[n, m]

                miu_pixel_distributions[m, j] = total * 1.0 / N_m

        return miu_pixel_distributions


    def maximize_pi(self,n_images, z_matrix, x_vectors, k_clusters, pi_proportion_of_images):
        for m in xrange(self.k_clusters):
            N_m = self.N_m_effective_number_of_images(z_matrix, m)
            pi_proportion_of_images[m] = N_m * 1.0 / n_images

        return pi_proportion_of_images


    def find_labels(self,z_matrix):
        labels = []

        for i in xrange(z_matrix.shape[0]):
            labels.append(argmax(z_matrix[i]))
        return labels


    def get_keys_with_more_than_one_el(self,dic):
        keys = []
        for key in dic.keys():
            if len(dic[key])>1:
                keys.append(key)
        return keys

    def fill_the_gaps(self,final_dict, k):

        for i in xrange(k):
            if i not in final_dict.keys():
                if len(final_dict.keys()) == 1:
                    key = final_dict.keys()[0]
                else:
                    possible_keys = self.get_keys_with_more_than_one_el(final_dict)
                    key = random.choice(possible_keys)

                final_dict[i] = [final_dict[key].pop(0)]

        return final_dict

    def do_clustering(self):

        x_vectors = []
        for key in self.unique_rows.keys():
            rep = self.unique_rows[key]
            x_vector = [str(el) for el in key.split(",")]
            for i in xrange(rep):
                x_vectors.append(x_vector)

        self.x_vectors = array(x_vectors)
        self.n_images = sum(self.unique_rows[key] for key in self.unique_rows.keys())
        self.d_pixels = len(self.unique_rows.keys()[0].split(","))



        # true_labels = []
        # for key in self.unique_rows.keys():
        #     label = get_label_of_cluster(vector=key, full_dict = self.full_data_dict)
        #     for _ in range(self.unique_rows[key]):
        #         true_labels.append(label)

        silhoutes_scores =[]
        rands = []
        import timeit
        start = timeit.default_timer()
        for i in range(self.number_of_iterations):
            miu_pixel_distributions = ones((self.k_clusters, self.d_pixels))

            for k in xrange(self.k_clusters):
                normalize = 0

                for i in xrange(self.d_pixels):
                    miu_pixel_distributions[k, i] *= (random_sample() * 0.5) + 0.25
                    normalize += miu_pixel_distributions[k, i]

                for i in xrange(self.d_pixels):
                    miu_pixel_distributions[k, i] /= normalize

            z_matrix = zeros((self.n_images, self.k_clusters), dtype=float64)

            z_matrix[:, 0] = ones(self.n_images)
            for i in range(self.n_images):
                shuffle(z_matrix[i])

            init_labels = self.find_labels(z_matrix)
            init_dic = {}

            for i in xrange(self.n_images):
                try:
                    init_dic[init_labels[i]].append(x_vectors[i])
                except:
                    init_dic[init_labels[i]] = []
                    init_dic[init_labels[i]].append(x_vectors[i])

            pi_proportion_of_images = sum(z_matrix, axis=0) * 1.0 / self.n_images

            for i in xrange(self.number_of_iterations):
                # Expectation step
                z_matrix = self.expectation_step(z_matrix, pi_proportion_of_images, miu_pixel_distributions)

                # Maximization step
                miu_pixel_distributions = self.maximize_miu(miu_pixel_distributions, z_matrix, x_vectors)
                pi_proportion_of_images = self.maximize_pi(self.n_images, z_matrix, x_vectors, self.k_clusters,
                                                      pi_proportion_of_images)

            final_dict = {}
            # print true_labels
            labels = self.find_labels(z_matrix)
            # print labels
            predicted_labels = labels

            for i in xrange(self.n_images):
                try:
                    final_dict[labels[i]].append(",".join(str(el) for el in self.x_vectors[i]))
                except:
                    final_dict[labels[i]] = [",".join(str(el) for el in self.x_vectors[i])]
            predicted_labels = asarray(predicted_labels)
            if len(np.unique(predicted_labels)) == 1:
                continue
            sil = silhouette_score(self.distance_matrix,predicted_labels, metric="precomputed")
            silhoutes_scores.append(sil)


            # rand = adjusted_rand_score(true_labels, predicted_labels)
            # rands.append(rand)

            # stop = timeit.default_timer()
            return final_dict
            # print "----SUMMARY---------- with iterations:", self.number_of_iterations
            # print "True Silhoutte score", silhouette_score(self.distance_matrix, true_labels, metric="precomputed")
            # print "Average silhoute score", np.mean(np.nan_to_num(np.asarray(silhoutes_scores)))
            # print "Average Adjusted rand index", np.mean(np.asarray(rands))
            # print "Time / sample", (stop - start) / float(n_times)

if __name__ == '__main__':
    set_printoptions(threshold=nan)
    # directory = "/home/laurynas/workspace/individual_project/data/sottoriva/CT_IRX2P_L2.csv"

    # unique_rows, full_or_data = process_file(directory)

    sample_name = "analysis_10_10_0.01_100"
    dir = "/home/laurynas/workspace/individual_project/simulated_data/"+sample_name+".txt"
    # ------------------------------------------------------------------------------------------------
    k_clusters = 10
    lam = 5
    loop_times = 1
    unique_rows, full_data_dict, full_info = read_simulated_data_file(dir)

    bmm = BMM(unique_rows,full_data_dict,full_info,number_of_iterations=20)

    for j in range(1,loop_times+1):
        bmm.do_clustering()


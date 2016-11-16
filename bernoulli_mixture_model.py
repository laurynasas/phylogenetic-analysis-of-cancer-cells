from numpy import *
from numpy.random import *

from k_medoid_clustering import *


def count_occurences(list):
    occ = {}
    for el in list:
        try:
            occ[el] += 1
        except:
            occ[el] = 1

    total = sum(occ[key] for key in occ.keys())

    for key in occ.keys():
        occ[key] = occ[key] * 1.0 / total
    print occ
    return occ


def bernoulli_distro(p, x):
    return power(p, x) * power((1 - p), (1 - x))


def prob_x_generated_by_cluster_k(x_vector, miu_pixel_distributions, k_cluster):
    product = 1

    for i in xrange(len(x_vector)):
        product *= bernoulli_distro(miu_pixel_distributions[k_cluster][i], x_vector[i])
    return product


def likelihood_of_single_training_image(x_vector, pi_proportion_of_images, miu_pixel_distributions):
    k_clusters = len(pi_proportion_of_images)

    likelihood = 0

    for k_cluster in xrange(k_clusters):
        likelihood += pi_proportion_of_images[k_cluster] * prob_x_generated_by_cluster_k(x_vector,
                                                                                         miu_pixel_distributions,
                                                                                         k_cluster)

    return likelihood


def likelihood_complete_data(pi_proportion_of_images, miu_pixel_distributions, x_vectors, z_matrix):
    n_images, d_pixels = x_vectors.shape

    k_clusters = len(pi_proportion_of_images)

    likelihood = 1
    for n in xrange(n_images):
        for k in xrange(k_clusters):
            likelihood *= power(
                pi_proportion_of_images[k] * prob_x_generated_by_cluster_k(x_vectors[n], miu_pixel_distributions, k),
                z_matrix[n, k])
    return likelihood


def log_likelihood_complete_data(likelihood):
    return log(likelihood)


def expectation_step(z_matrix, pi_proportion_of_images, miu_pixel_distributions):
    n_images, k_clusters = z_matrix.shape

    for n in xrange(n_images):
        x_vector = x_vectors[n]
        for k in xrange(k_clusters):
            nom = pi_proportion_of_images[k] * prob_x_generated_by_cluster_k(x_vector, miu_pixel_distributions, k)
            den = 0
            for m in xrange(k_clusters):
                den += pi_proportion_of_images[m] * prob_x_generated_by_cluster_k(x_vector, miu_pixel_distributions, m)

            z_matrix[n, k] = nom * 1.0 / den

    return z_matrix


def N_m_effective_number_of_images(z_matrix, k_cluster):
    total = 0
    for n in xrange(n_images):
        total += z_matrix[n, k_cluster]
    return total


def maximize_miu(miu_pixel_distributions, z_matrix, x_vectors):

    for m in xrange(k_clusters):
        for j in xrange(d_pixels):
            N_m = N_m_effective_number_of_images(n_images, z_matrix, m)

            total = sum(x_vectors[n, j] * z_matrix[n, m] for n in xrange(n_images))

            miu_pixel_distributions[m, j] = total * 1.0 / N_m

    return miu_pixel_distributions


def maximize_pi(n_images, z_matrix, x_vectors, k_clusters, pi_proportion_of_images):
    for m in xrange(k_clusters):
        N_m = N_m_effective_number_of_images(n_images, z_matrix, m)
        pi_proportion_of_images[m] = N_m * 1.0 / n_images

    return pi_proportion_of_images


if __name__ == '__main__':
    set_printoptions(threshold=nan)
    directory = "./data/sample_data.csv"
    unique_rows = process_file(directory)
    global n_images, d_pixels, k_clusters, x_vectors

    x_vectors = []
    for key in unique_rows.keys():
        rep = unique_rows[key]
        x_vector = [int(el) for el in key.split(",")]
        for i in xrange(rep):
            x_vectors.append(x_vector)

    x_vectors = array(x_vectors)
    n_images = sum(unique_rows[key] for key in unique_rows.keys())
    d_pixels = len(unique_rows.keys()[0].split(","))
    k_clusters = 12

    miu_pixel_distributions = zeros((d_pixels, k_clusters))

    lam = 5
    # pi_proportion_of_images = dirichlet(ones(k_clusters), size=1)

    z_matrix = zeros((k_clusters, n_images), dtype=int64)

    z_matrix[0] = ones(n_images)

    for i in range(n_images):
        shuffle(z_matrix[:, i])
    print sum(z_matrix, axis=1)
    pi_proportion_of_images = sum(z_matrix, axis=1) * 1.0 / n_images

    # print z_matrix
    # pi_proportion_of_images = pi_proportion_of_images*1.0 / k_clusters
    # occ_dict = count_occurences(pi_proportion_of_images)

    # print pi_proportion_of_images
    # print unique_rows
    # M, C = do_k_medoid_clustering(distance_matrix, 10)
    # print type(distance_matrix.shape[0])
    # labels = np.zeros(distance_matrix.shape[0])
    # print C
    # print unique_rows
    # for key in zip(C.keys()):
    #     value_list = C[key[0]]
    #     for value in value_list:
    #         labels[value] = key[0]
    #
    # print silhouette_score_slow(distance_matrix, labels)

import itertools

from numpy.random import *

from k_medoid_clustering import *
from rand_index import *

RESTRICTED_CHARS = ["NA\n", "NA"]


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
    return occ


def bernoulli_distro(p, x):
    return (p ** x) * ((1 - p) ** (1 - x))


def prob_x_generated_by_cluster_k(x_vector, miu_pixel_distributions, k_cluster):
    product = 1
    for i in xrange(len(x_vector)):
        if x_vector[i] not in RESTRICTED_CHARS:
            pwr = int(x_vector[i])
            product *= bernoulli_distro(miu_pixel_distributions[k_cluster][i], pwr)
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
    for n in xrange(n_images):
        x_vector = x_vectors[n]
        for k in xrange(k_clusters):
            nom = pi_proportion_of_images[k] * prob_x_generated_by_cluster_k(x_vector, miu_pixel_distributions, k)
            den = 0
            for m in xrange(k_clusters):
                den += pi_proportion_of_images[m] * prob_x_generated_by_cluster_k(x_vector, miu_pixel_distributions, m)

            nom = nan_to_num(nom)
            den = nan_to_num(den)

            # if den == 0.0:
            #     den = 1.0
            # print nom, den
            z_matrix[n, k] = nom * 1.0 / den
    # print z_matrix, "\n------^z_new---\n"
    return z_matrix


def N_m_effective_number_of_images(z_matrix, k_cluster):
    total = 0
    for n in xrange(n_images):
        total += z_matrix[n, k_cluster]
    return total


def maximize_miu(miu_pixel_distributions, z_matrix, x_vectors):
    for m in xrange(k_clusters):
        for j in xrange(d_pixels):
            N_m = N_m_effective_number_of_images(z_matrix, m)

            total = 0
            for n in xrange(n_images):
                if x_vectors[n, j] not in RESTRICTED_CHARS:
                    total += int(x_vectors[n, j]) * z_matrix[n, m]

            miu_pixel_distributions[m, j] = total * 1.0 / N_m

    return miu_pixel_distributions


def maximize_pi(n_images, z_matrix, x_vectors, k_clusters, pi_proportion_of_images):
    for m in xrange(k_clusters):
        N_m = N_m_effective_number_of_images(z_matrix, m)
        pi_proportion_of_images[m] = N_m * 1.0 / n_images

    return pi_proportion_of_images


def find_labels(z_matrix):
    labels = []

    for i in xrange(z_matrix.shape[0]):
        labels.append(argmax(z_matrix[i]))
    return labels


def get_keys_with_more_than_one_el(dic):
    keys = []
    for key in dic.keys():
        if len(dic[key])>1:
            keys.append(key)
    return keys

def fill_the_gaps(final_dict, k):
    # print "there are keys: ", final_dict.keys()
    for i in xrange(k):
        # print sum([len(final_dict[key]) for key in final_dict.keys()])
        if i not in final_dict.keys():
            # print "this is not in keys needs to be filled", i
            if len(final_dict.keys()) == 1:
                key = final_dict.keys()[0]
            else:
                possible_keys = get_keys_with_more_than_one_el(final_dict)
                key = random.choice(possible_keys)

            # print "picked key to take elemnts from ", key

            final_dict[i] = [final_dict[key].pop(0)]
    # print "populated",final_dict
    return final_dict


if __name__ == '__main__':
    set_printoptions(threshold=nan)
    # directory = "/home/laurynas/workspace/individual_project/data/sottoriva/CT_IRX2P_L2.csv"

    # unique_rows, full_or_data = process_file(directory)

    sample_name = "analysis_20_20_0.01_1000"
    dir = "/home/laurynas/workspace/individual_project/simulated_data/"+sample_name+".txt"
    n_times = 10
    # ------------------------------------------------------------------------------------------------
    k_clusters = 20
    lam = 5
    iterations = 25
    for j in range(1,iterations+1):
        iter_count = j

        # ------------------------------------------------------------------------------------------------

        # directory = "/home/laurynas/workspace/individual_project/data/hou/snv_cut.csv"
        # unique_rows, full_or_data = process_single_cell_data_file(directory)
        unique_rows, full_data_dict, full_info = read_simulated_data_file(dir)
        distance_matrix = np.matrix(find_distance_matrix(unique_rows))

        true_labels = []
        for key in unique_rows.keys():
            label = get_label_of_cluster(vector=key, full_dict = full_data_dict)
            for _ in range(unique_rows[key]):
                true_labels.append(label)

        # print true_labels
        global n_images, d_pixels, k_clusters, x_vectors

        x_vectors = []
        for key in unique_rows.keys():
            rep = unique_rows[key]
            x_vector = [str(el) for el in key.split(",")]
            for i in xrange(rep):
                x_vectors.append(x_vector)

        # full_or_data = {1:[[1, 0, 1, 1], [1, 1, 1, 1]], 2: [[0, 0, 0, 1], [0, 0, 0, 0]], 3: [[0, 0, 1, 1]]}
        # x_vectors = [[1, 0, 1, 1], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1], [1, 1, 1, 1]]
        x_vectors = array(x_vectors)

        n_images = sum(unique_rows[key] for key in unique_rows.keys())
        # n_images = len(x_vectors)

        d_pixels = len(unique_rows.keys()[0].split(","))
        # d_pixels = len(x_vectors[0])


        silhoutes_scores =[]
        rands = []
        import timeit
        start = timeit.default_timer()
        for i in range(n_times):
            # print i
            miu_pixel_distributions = ones((k_clusters, d_pixels))

            for k in xrange(k_clusters):
                normalize = 0

                for i in xrange(d_pixels):
                    miu_pixel_distributions[k, i] *= (random_sample() * 0.5) + 0.25
                    normalize += miu_pixel_distributions[k, i]

                for i in xrange(d_pixels):
                    miu_pixel_distributions[k, i] /= normalize

            z_matrix = zeros((n_images, k_clusters), dtype=float64)

            z_matrix[:, 0] = ones(n_images)
            for i in range(n_images):
                shuffle(z_matrix[i])

            init_labels = find_labels(z_matrix)
            init_dic = {}

            for i in xrange(n_images):
                try:
                    init_dic[init_labels[i]].append(x_vectors[i])
                except:
                    init_dic[init_labels[i]] = []
                    init_dic[init_labels[i]].append(x_vectors[i])


            # print "\n-------------------------------\n"

            pi_proportion_of_images = sum(z_matrix, axis=0) * 1.0 / n_images

            for i in xrange(iter_count):
                # Expectation step
                z_matrix = expectation_step(z_matrix, pi_proportion_of_images, miu_pixel_distributions)

                # Maximization step
                miu_pixel_distributions = maximize_miu(miu_pixel_distributions, z_matrix, x_vectors)
                pi_proportion_of_images = maximize_pi(n_images, z_matrix, x_vectors, k_clusters, pi_proportion_of_images)

            final_dict = {}
            # print true_labels
            labels = find_labels(z_matrix)
            # print labels
            predicted_labels = labels
            for i in xrange(n_images):
                try:
                    final_dict[labels[i]].append(x_vectors[i].tolist())
                except:
                    final_dict[labels[i]] = [x_vectors[i].tolist()]

            predicted_labels = asarray(predicted_labels)
            if len(np.unique(predicted_labels)) == 1:
                continue
            sil = silhouette_score(distance_matrix,predicted_labels, metric="precomputed")
            silhoutes_scores.append(sil)


            rand = adjusted_rand_score(true_labels, predicted_labels)
            rands.append(rand)

            # for key in final_dict.keys():
            #     # for el in final_dict[key]:
            #         # pass
            #     print "Cluster", key," and el | ", final_dict[key]
            #         # pass

            # print "No of clusters before filling", len(final_dict.keys())
            # full_dict = fill_the_gaps(final_dict, k_clusters)
            # print "No of clusters after filling", len(full_dict.keys())


            # def product_gen(n):
            #     for r in itertools.count(1):
            #         for x in itertools.product(n, repeat=r):
            #             yield "".join(x)



            # from parsimony_tree_builder import get_genotypes_from_clusters
            # f = open('./data/genotype.phy', 'w')
            # print full_dict
            # genotypes = get_genotypes_from_clusters(full_dict, d_pixels)
            # print '--->', full_dict.keys(),len(genotypes)
            # str_genotype = ''
            # for label, gen in zip(product_gen(string.ascii_uppercase), genotypes):
            #     str1 = ''.join([str(el) for el in gen])
            #     str_genotype += str(label).ljust(10)+str1+'\n'
            # f.write(str(len(genotypes))+' '+str(d_pixels)+'\n')
            # f.write(str_genotype)
            # f.close()


            # print "Rand index: ", rand_index(full_or_data, full_dict)
            # labeling = find_common_labelling(full_or_data, final_dict)
            # calculate_tp_fp(labeling, full_or_data, full_dict)

            # print miu_pixel_distributions
            # print "after",z_matrix
            # print sum(z_matrix, axis=1)
            # print pi_proportion_of_images

            #
            # from processing_sample_data import write_to_file_and_save_image
            #
            # write_to_file_and_save_image("bernoulli_mixture_model", "snv_cut_biger", final_dict, k_clusters)

            # predicted_labels = asarray(predicted_labels)
            # print "silh: ", silhouette_score_slow(distance_matrix, asarray(true_labels), metric="precomputed")
            # print "theor silh", silhouette_score(distance_matrix, predicted_labels, metric="precomputed")
            #
            # print "Adjusted rand: ", adjusted_rand_score(true_labels, predicted_labels)
            # print "Will be saving the image", plot_2D_similarity_matrix(distance_matrix, "K_medoids", "analysis_5_5_0.05_20",
            #                                                             no_cl=str(k_clusters), data_type="simulated")
        stop = timeit.default_timer()
        print "----SUMMARY---------- with iterations:",iter_count
        print "True Silhoutte score", silhouette_score(distance_matrix,true_labels, metric="precomputed")
        print "Average silhoute score", np.mean(np.nan_to_num(np.asarray(silhoutes_scores)))
        # print "Average Adjusted rand index", np.mean(np.asarray(rands))
        print "Time / sample", (stop - start)/float(n_times)
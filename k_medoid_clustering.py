from sklearn.cluster import k_means_
from sklearn.metrics import adjusted_rand_score, silhouette_score

from parsimony_tree_builder import get_genotypes_from_clusters
from processing_sample_data import *
from silhouette_score_implementation import *
from single_linkage_clustering import get_label_of_cluster
from single_linkage_clustering import read_simulated_data_file


def do_k_medoid_clustering(D, k, tmax=100):
    # determine dimensions of distance matrix D
    m, n = D.shape

    # randomly initialize an array of k medoid indices
    M = np.sort(np.random.choice(n, k))

    # create a copy of the array of medoid indices
    Mnew = np.copy(M)

    # initialize a dictionary to represent clusters
    C = {}
    for t in xrange(tmax):
        # determine clusters, i. e. arrays of data indices
        J = np.argmin(D[:, M], axis=1)
        for kappa in range(k):
            C[kappa] = np.where(J == kappa)[0]
        # update cluster medoids
        for kappa in range(k):
            # print D[np.ix_(C[kappa],C[kappa])]

            J = np.mean(D[np.ix_(C[kappa], C[kappa])], axis=1)
            if len(J) != 0:
                j = np.argmin(J)
                Mnew[kappa] = C[kappa][j]
        np.sort(Mnew)
        # check for convergence
        if np.array_equal(M, Mnew):
            break
        M = np.copy(Mnew)
    else:
        # final update of cluster memberships
        J = np.argmin(D[:, M], axis=1)
        for kappa in range(k):
            C[kappa] = np.where(J == kappa)[0]

    # return results
    return M, C


def compare_genotypes(true_genotype, predicted_genotype):
    count = 0
    for gen in predicted_genotype:
        if gen in true_genotype:
            count += 1

    return count / float(len(true_genotype))


def get_true_genotypes(dir):
    simulated_data_file = open(dir, 'rw+')
    data = simulated_data_file.readlines()
    true_genotypes = []
    for line in data:
        true_genotypes.append(map(int, line.split(',')))

    return true_genotypes


if __name__ == '__main__':
    sample_name = "analysis_for_genotypes_20_20_0.1_1000"
    dir = "/home/laurynas/workspace/individual_project/simulated_data/" + sample_name + ".txt"
    n_times = 1
    no_cl = 20
    vector_size = 20
    unique_rows, full_data_dict, full_info = read_simulated_data_file(dir)

    sample_name = "analysis_genotypes_20_20_0.1_1000"
    dir = "/home/laurynas/workspace/individual_project/simulated_data/" + sample_name + ".txt"

    true_genotypes = get_true_genotypes(dir)

    distance_matrix = np.matrix(find_distance_matrix(unique_rows))
    print distance_matrix
    unique_strings = unique_rows.keys()
    # print unique_rows

    true_labels = []
    for key in unique_rows.keys():
        label = get_label_of_cluster(vector=key, full_dict=full_data_dict)
        for _ in range(unique_rows[key]):
            true_labels.append(label)
    print true_labels

    # print "->",full_data_dict
    predicted_labels = [0 for _ in range(len(full_info))]
    # print predicted_labels

    silhoutes_scores = []
    rands = []
    import timeit

    start = timeit.default_timer()
    for i in range(n_times):
        print i

        M, C = do_k_medoid_clustering(distance_matrix, no_cl)
        k = k_means_.KMeans(n_clusters=no_cl, random_state=0).fit(distance_matrix)
        full_dict = {}
        for label, vector in zip(k.labels_, unique_rows.keys()):
            if label not in full_dict.keys():
                full_dict[label] = [vector]
            else:
                full_dict[label] += [vector]

        # counter = 0
        # for key,value in C.items():
        #     for el in value:
        #         # print el
        #         predicted_labels[counter] = key
        #         counter +=1

        # predicted_labels = asarray(predicted_labels)
        predicted_labels = k.labels_
        if len(np.unique(predicted_labels)) == 1:
            continue
        sil = silhouette_score(distance_matrix, predicted_labels, metric="precomputed")
        silhoutes_scores.append(sil)

        rand = adjusted_rand_score(true_labels, predicted_labels)
        rands.append(rand)
        # print full_dict
        predicted = get_genotypes_from_clusters(full_dict, vector_size)

        print "True genotypes: ", true_genotypes
        print "Predicted genotypes: ",predicted
        print "Compared genotypes correct/total_true: ", compare_genotypes(true_genotypes, predicted)
        # print "silh: ",silhouette_score_slow(distance_matrix, asarray(true_labels), metric="precomputed")
        # print "theor silh", silhouette_score(distance_matrix, predicted_labels, metric="precomputed")

        # print "Adjusted rand: ", adjusted_rand_score(true_labels, predicted_labels)
        # print "Will be saving the image", plot_2D_similarity_matrix(distance_matrix,"K_medoids","analysis_5_5_0.05_20",no_cl=str(no_cl),data_type="simulated")
    stop = timeit.default_timer()
    print "----SUMMARY----------"
    print "True Silhoutte score", silhouette_score(distance_matrix, true_labels, metric="precomputed")
    print "Average silhoute score", np.mean(np.nan_to_num(np.asarray(silhoutes_scores)))
    print "Average Adjusted rand index", np.mean(np.asarray(rands))
    print "Time / sample", (stop - start) / float(n_times)

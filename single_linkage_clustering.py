import numpy as np
from scipy.cluster.hierarchy import fcluster,linkage
import processing_sample_data as pr
from silhouette_score_implementation import silhouette_score_slow
from rand_index import rand_index
from sklearn.metrics import adjusted_rand_score,silhouette_score
import parsimony_tree_builder

class DataNode:
    def __init__(self, vector, cluster_label):
        self.vector = vector
        self.cluster_label = cluster_label

def read_simulated_data_file(dir):
    simulated_data_file = open(dir, 'rw+')
    unique_rows = {}
    full_data_dict = {}
    full_info =[]
    for line in simulated_data_file:
        # line = line[:-1]
        # label = int(line.split(",")[-1])
        line = line.split(" | ")
        cluster_label = int(line[1])
        vector = line[0][1:-1]
        full_info.append(DataNode(vector,cluster_label))

        if full_data_dict.get(cluster_label):
            full_data_dict[cluster_label] += [map(int, vector.split(","))]
        else:
            full_data_dict[cluster_label] = [map(int, vector.split(","))]
        if vector in unique_rows.keys():
            unique_rows[vector] += 1
        else:
            unique_rows[vector] = 1

    return unique_rows, full_data_dict, full_info

def find_min_corner(distance_matrix):
    # print distance_matrix
    minimal = max(distance_matrix[0][1:])
    location = (-1, -1)
    for i in range(len(distance_matrix[0])):
        for j in range(i, len(distance_matrix[0])):
            if distance_matrix[i][j] != 0 and distance_matrix[i][j] < minimal:
                minimal = distance_matrix[i][j]
                location = (i, j)
    print minimal, location
    return location


def merge(distance_matrix, location, clusters):
    # print len(distance_matrix[0])
    distance_matrix = distance_matrix.tolist()
    new_row = [min(distance_matrix[location[0]][i], distance_matrix[location[1]][i]) for i in
               range(len(distance_matrix[0]))]
    del new_row[location[1]]
    print "new row: ", new_row
    distance_matrix = np.matrix(distance_matrix)
    print "Cluster before failing", clusters
    clusters[location[0]] += [clusters[location[1]]]
    del clusters[location[1]]
    print "Cluster after adding new stuff", clusters
    distance_matrix = np.delete(distance_matrix, (location[0]), axis=0)
    distance_matrix = np.delete(distance_matrix, (location[1] - 1), axis=0)

    distance_matrix = np.delete(distance_matrix, (location[0]), axis=1)
    distance_matrix = np.delete(distance_matrix, (location[1] - 1), axis=1)
    distance_matrix = np.asarray(distance_matrix)
    print "Before adding new lines\n", distance_matrix
    distance_matrix = np.insert(distance_matrix, location[0], new_row[:location[0]] + new_row[location[0] + 1:], axis=0)
    print "Before adding new lines\n", distance_matrix

    distance_matrix = np.insert(distance_matrix, location[0], new_row, axis=1)
    distance_matrix = np.matrix(distance_matrix)
    print "After added new liens\n", distance_matrix
    print "SIZE after merging", distance_matrix.size

    return distance_matrix, clusters


# print new_row
# print distance_matrix
# delete_row_and_column(distance_matrix, location)
# add_new_row_and_column(distance_matrix, new_row, location)

def single_linkage_clustering(distance_matrix):
    clusters = []
    for i in range(distance_matrix[0].size):
        clusters += [['V' + str(i)]]
    # print clusters
    location = find_min_corner(distance_matrix.tolist())
    while location != (-1, -1):
        location = find_min_corner(distance_matrix.tolist())
        if location != (-1, -1):
            (distance_matrix, clusters) = merge(distance_matrix, location, clusters)
    return clusters
def get_label_of_cluster(vector, full_dict):
    vector = map(int, vector.split(","))
    for key in full_dict.keys():
        # print vector?S, full_dict[key]
        if vector in full_dict[key]:
            return key

if __name__ == "__main__":
    sample_name = "populated_true_genotypes_10_10_0.01_100"
    dir = "/home/laurynas/workspace/individual_project/simulated_data/"+sample_name+".txt"
    n_times = 1

    unique_rows, full_data_dict, full_info = read_simulated_data_file(dir)
    distance_matrix = np.matrix(pr.find_distance_matrix(unique_rows))
    # print distance_matrix
    true_labels = []
    for key in unique_rows.keys():
        label = get_label_of_cluster(vector=key, full_dict = full_data_dict)
        for _ in range(unique_rows[key]):
            true_labels.append(label)

    silhoutes_scores =[]
    rands = []
    import timeit
    start = timeit.default_timer()
    for i in range(n_times):
        print i
        labels = fcluster(linkage(distance_matrix, method='complete'), t=15, criterion='maxclust')
        # print labels-1, unique_rows.keys()
        # print unique_rows.keys()
        distribution= {}
        for key, label in zip(unique_rows.keys(),labels-1):
            if label in distribution.keys():
                distribution[label] += [key]
            else:
                distribution[label] = [key]
        new_data = {}

        for cluster_label in labels-1:
            if new_data.get(cluster_label):
                new_data[cluster_label] += [map(int, unique_rows.keys()[cluster_label].split(","))]
            else:
                new_data[cluster_label] = [map(int, unique_rows.keys()[cluster_label].split(","))]
        new_data_formatetd = {}
        for key,value in new_data.items():
            new_value=[','.join(map(str,el)) for el in value]
            new_data_formatetd[key] = new_value
        # print new_data_formatetd
        ready_genotype = parsimony_tree_builder.get_genotypes_from_clusters(new_data_formatetd,10)
        print ready_genotype
        sample_name = "predicted_SLC_ready_genotypes_10_10_0.01_100"
        dir = "/home/laurynas/workspace/individual_project/simulated_data/" + sample_name + ".phy"
        parsimony_tree_builder.write_genotype_to_file(dir, ready_genotype)

        # for key in distribution:
        #     for el in distribution[key]:
        #         for _ in xrange(unique_rows[str(el)]):
        #
        #             print "[",str(el),"]", " | ", key

        vectors =[distribution[key] for key in distribution.keys()]
        # print distance_matrix
        # distance_matrix = np.ones((len(labels),len(labels)))
        sil = silhouette_score(distance_matrix,labels-1, metric="precomputed")
        silhoutes_scores.append(sil)
        # print "predicted silhoutte score", sil
        # print rand_index(full_data_dict, new_data)
        labels_predicted = labels-1
        # print true_labels,labels_predicted
        rand = adjusted_rand_score(true_labels, labels_predicted)
        rands.append(rand)
        # print "Adjusted rand: ", rand
        # print "Will be saving the image", plot_2D_similarity_matrix(distance_matrix,"SLC","gdrive_"+sample_name,data_type="simulated")
    stop = timeit.default_timer()
    print "----SUMMARY----------"
    print "True Silhoutte score", silhouette_score(distance_matrix,true_labels, metric="precomputed")
    print "Average silhoute score", np.mean(np.asarray(silhoutes_scores))
    print "Average Adjusted rand index", np.mean(np.asarray(rands))
    print "Time / sample", (stop - start)/float(n_times)

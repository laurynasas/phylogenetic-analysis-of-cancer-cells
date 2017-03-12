from numpy import *
from numpy.random import choice
import datetime

def generate_vector(p, len_vectors, parent_vector):
    new_vector = parent_vector[:]
    for i in xrange(len_vectors):
        if parent_vector[i] == 0:
            choice_one = choice(array([0, 1]), p=[1 - p, p])
            # print choice_one
            new_vector[i] = choice_one
    return new_vector


def generate_original_clusters(k_clusters, len_vectors):
    tree = [[]]
    p_increment = 1 / float(k_clusters)
    root_vector = [0] * len_vectors
    tree.append(root_vector)
    p = p_increment
    for i in xrange(2, k_clusters + 1):
        if i % 2 == 1:
            parent_vector = tree[(i - 1) / 2]
        else:
            parent_vector = tree[i / 2]
        child = generate_vector(p, len_vectors, parent_vector)
        while child in tree:
            child = generate_vector(p, len_vectors, parent_vector)

        tree.append(child)
        p += p_increment
    return tree[1:]


def get_cluster_sizes(k_clusters, dataset_size, low=1):
    if dataset_size < k_clusters:
        return None
    high = dataset_size
    # print high
    possible_numbers = range(low, high)
    # print possible_numbers
    sizes = []
    for _ in xrange(k_clusters - 1):
        if possible_numbers:
            picked_number = choice(possible_numbers)
            sizes.append(picked_number)
            possible_numbers.pop(possible_numbers.index(picked_number))
        else:
            break

    sizes = sorted(sizes)

    cluster_sizes = [sizes[0]]
    for i in xrange(k_clusters - 2):
        cluster_sizes.append(sizes[i + 1] - sizes[i])
    cluster_sizes.append(dataset_size - sizes[-1])
    return cluster_sizes


def flip(err, el):
    if el == 1:
        return choice(array([0, 1]), p=[err, 1 - err])
    if el == 0:
        return choice(array([1, 0]), p=[err, 1 - err])
    return None


def get_vector_with_error(vector, error):
    new_vector = []
    for el in vector:
        new_vector.append(flip(error, el))
    return new_vector


def populate_cluster_with_errorous_data(orig_data, cluster_sizes, error):
    labeled_data = {}
    for index, ideal_el in enumerate(orig_data):
        labeled_data[index] = []
        for j in xrange(cluster_sizes[index]):
            labeled_data[index].append(get_vector_with_error(ideal_el, error))

    return labeled_data


k_clusters = 10
len_vectors = 10
error_perct = 0.01
dataset_size = 100

original_clusters = generate_original_clusters(k_clusters, len_vectors)

cluster_sizes = get_cluster_sizes(k_clusters, dataset_size)
populated =  populate_cluster_with_errorous_data(original_clusters, cluster_sizes, error_perct)

print original_clusters
target = open("./simulated_data/true_genotypes_"+str(len_vectors)+"_"+str(k_clusters)+"_"+str(error_perct)+"_"+str(dataset_size)+".txt", 'w+')
for genotype in original_clusters:
    target.write(",".join(map(str,genotype)) +"\n")

target.close()


target = open("./simulated_data/populated_true_genotypes_"+str(len_vectors)+"_"+str(k_clusters)+"_"+str(error_perct)+"_"+str(dataset_size)+".txt", 'w+')
for key in populated.keys():
    for el in populated[key]:
        target.write(str(el) + ' | ' +str(key)+"\n")

target.close()
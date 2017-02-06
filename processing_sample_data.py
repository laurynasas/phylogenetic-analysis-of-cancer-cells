from scipy.cluster.hierarchy import fcluster, linkage

from rand_index import *
from single_linkage_clustering import *


def process_file(dir):
    sample_data_file = open(dir, 'rw+')
    unique_rows = {}
    full_data_dict = {}
    sample_data_file.readline()

    for line in sample_data_file:
        line = line[:-1]
        label = int(line.split(",")[-1])
        line = ','.join([str(x) for x in (line.split(",")[:-1])])
        # print "label",label
        # label = int(line.split(",")[-1])
        if full_data_dict.get(label):
            full_data_dict[label] += [map(int, line.split(","))]
        else:
            full_data_dict[label] = [map(int, line.split(","))]
        if line in unique_rows.keys():
            unique_rows[line] += 1
        else:
            unique_rows[line] = 1

    return unique_rows, full_data_dict

def process_single_cell_data_file(dir):
    sample_data_file = open(dir, 'rw+')
    unique_rows = {}
    full_data_dict = {}
    sample_data_file.readline()

    for line in sample_data_file:
        # line = line[:-1]
        label = line.split('"')[1]
        line = (line.split('"')[2])[1:]
        # print "label",label
        # print line
        # label = int(line.split(",")[-1])
        if full_data_dict.get(label):
            full_data_dict[label] += [map(str, line.split(","))]
        else:
            full_data_dict[label] = [map(str, line.split(","))]
        if line in unique_rows.keys():
            unique_rows[line] += 1
        else:
            unique_rows[line] = 1

    return unique_rows, full_data_dict

def diff_char_string(a, b):
    distance = 0
    for i, j in zip(a, b):
        if i != j:
            distance += 1
    return distance


def find_distance_matrix(unique_rows):
    number_rows = len(unique_rows.keys())
    unique_strings = unique_rows.keys()

    distance_matrix = []
    for i in range(number_rows):
        new_row = []
        for j in range(number_rows):
            new_row += [diff_char_string(unique_strings[i], unique_strings[j])]
        distance_matrix += [new_row]

    return distance_matrix


if __name__ == '__main__':
    # ----- NORMAL PROCESSING -------------------------------------------------------------------------------
    # dir = "./data/data2.txt"
    # unique_rows, full_data_dict = process_file(dir)
    # # print unique_rows
    # distance_matrix = np.matrix(find_distance_matrix(unique_rows))
    # # distance_matrix = np.matrix([[0,662,877,255,412,996],[662,0,295,468,268,400],[877,295,0,754,564,138],[255,468,754,0,219,869],[412,268,564,219,0,669],[996,400,138,869,669,0]])
    # # print "--->\n",distance_matrix
    # # print single_linkage_clustering(distance_matrix)
    #
    # labels = fcluster(linkage(distance_matrix, method='complete'), t=13, criterion='maxclust')
    # # print unique_rows
    # # print unique_rows
    # print labels
    # final_dict = {}
    # for index, label in enumerate(labels):
    #     print unique_rows.keys()[index], label
    #     if final_dict.get(label - 1):
    #         final_dict[label - 1] += [map(int, unique_rows.keys()[index].split(','))] * unique_rows.values()[index]
    #     else:
    #         final_dict[label - 1] = [map(int, unique_rows.keys()[index].split(','))] * unique_rows.values()[index]
    # print rand_index(full_data_dict, final_dict)
# ------------------------------------------------------------------------------------------------------------------

# -------------------------------------------SINGLE CELL PROCESSING-------------------------------------------------

    dir = "./data/hou/snv.csv"
    unique_rows, full_data_dict = process_single_cell_data_file(dir)
    print unique_rows

    # distance_matrix = np.matrix(find_distance_matrix(unique_rows))
    # # distance_matrix = np.matrix([[0,662,877,255,412,996],[662,0,295,468,268,400],[877,295,0,754,564,138],[255,468,754,0,219,869],[412,268,564,219,0,669],[996,400,138,869,669,0]])
    # # print "--->\n",distance_matrix
    # # print single_linkage_clustering(distance_matrix)
    #
    # labels = fcluster(linkage(distance_matrix, method='complete'), t=13, criterion='maxclust')
    # # print unique_rows
    # # print unique_rows
    # print labels
    # final_dict = {}
    # for index, label in enumerate(labels):
    #     print unique_rows.keys()[index], label
    #     if final_dict.get(label - 1):
    #         final_dict[label - 1] += [map(int, unique_rows.keys()[index].split(','))] * unique_rows.values()[index]
    #     else:
    #         final_dict[label - 1] = [map(int, unique_rows.keys()[index].split(','))] * unique_rows.values()[index]
    # print rand_index(full_data_dict, final_dict)


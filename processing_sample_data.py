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
    dir = "./data/sample_data.csv"
    unique_rows = process_file(dir)
    # print unique_rows
    distance_matrix = np.matrix(find_distance_matrix(unique_rows))
    # distance_matrix = np.matrix([[0,662,877,255,412,996],[662,0,295,468,268,400],[877,295,0,754,564,138],[255,468,754,0,219,869],[412,268,564,219,0,669],[996,400,138,869,669,0]])
    # print "--->\n",distance_matrix
    print single_linkage_clustering(distance_matrix)

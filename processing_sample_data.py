from rand_index import *
from datetime import datetime
from single_cell_processing import calc_distance_matrix_from_file, plot_2D_similarity_matrix

def write_to_file_and_save_image(method, dataset, final_dict, no_cl):
    dir = "./data/hou/clustered_data/" + str(method) + "_" + str(dataset) + str(datetime.now()) + ".txt"
    target = open(dir, 'w+')
    for key in final_dict.keys():
        for el in final_dict[key]:
            target.write(str(key) + ' | ' + ','.join(str(x) for x in el) + "\n")

    target.close()

    matrix = calc_distance_matrix_from_file(dir)
    plot_2D_similarity_matrix(matrix, method, dataset, no_cl)


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
    class DataNode:
        def __init__(self, vector, cluster_label):
            self.vector = vector
            self.cluster_label = cluster_label

    sample_data_file = open(dir, 'rw+')
    unique_rows = {}
    full_data_dict = {}
    full_info=[]
    labels = (sample_data_file.readline()).split(",")
    data_lines = [[] for _ in range(len(labels)-1)]
    for line in sample_data_file:
        i=0
        # line_no +=1
        # label = labels[line_no]
        # line = line[:-1]
        # label = line.split('"')[1]
        line = (line.split('"')[2][1:])
        # print line
        line = map(str, line.split(","))
        # print "line len", len(line), line
        # print "->",line
        full_info.append(DataNode(','.join(str(x) for x in line), None))
        for i in xrange(len(line)):
            data_lines[i].append(line[i])
        # print data_lines
        # print "label",label
        # print line
        # label = int(line.split(",")[-1])
        # if full_data_dict.get(label):
        #     full_data_dict[label] += [map(str, line.split(","))]
        # else:
        #     full_data_dict[label] = [map(str, line.split(","))]
        # if line in unique_rows.keys():
        #     unique_rows[line] += 1
        # else:
        #     unique_rows[line] = 1

    for label, column in zip(labels, data_lines):
        column = [x.replace("\n","") for x in column]
        full_data_dict[label] = column
        # print str(column)[1:-1]
        row = ','.join(str(x) for x in column)
        unique_rows[row] = 1
    # print unique_rows

    return unique_rows, full_data_dict, full_info


def diff_char_string(a, b):
    distance = 0
    for i, j in zip(a, b):
        if i != j:
            distance += 1
    return distance


def find_distance_matrix(unique_rows):
    number_rows = len(unique_rows.keys())
    unique_strings = []
    # unique_strings = unique_rows.keys()

    for key in unique_rows.keys():
        for _ in xrange(unique_rows[key]):
            unique_strings += [key]

    number_rows = len(unique_strings)

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

    # dir = "./data/hou/snv.csv"
    # unique_rows, full_data_dict = process_single_cell_data_file(dir)



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
    pass
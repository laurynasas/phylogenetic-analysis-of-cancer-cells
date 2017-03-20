'''
    Class was implemented to avoid circular imports and unify common methods in one class.
'''


class Helper:
    def __init__(self):
        pass

    '''
        Given some vector, returns its cluster label
    '''

    def get_label_of_cluster(self, vector, full_dict):
        vector = map(int, vector.split(","))
        for key in full_dict.keys():
            if vector in full_dict[key]:
                return key

    '''
        Gets Hammington distance between to binary vectors
    '''

    def diff_char_string(self, a, b):
        distance = 0
        for i, j in zip(a, b):
            if i != j:
                distance += 1
        return distance

    '''
        Return distance matrix based on the unique genotypes
    '''

    def find_distance_matrix(self, unique_rows):
        unique_strings = []

        for key in unique_rows.keys():
            for _ in xrange(unique_rows[key]):
                unique_strings += [key]

        number_rows = len(unique_strings)

        distance_matrix = []
        for i in range(number_rows):
            new_row = []
            for j in range(number_rows):
                new_row += [self.diff_char_string(unique_strings[i], unique_strings[j])]
            distance_matrix += [new_row]

        return distance_matrix

    '''
        Data node is auxiliary class to encapsulate all of the vector's information.
    '''

    class DataNode:
        def __init__(self, vector, cluster_label):
            self.vector = vector
            self.cluster_label = cluster_label

    '''
        Given raw simulated file directory it reads and processes file.
    '''

    def read_simulated_data_file(self, dir):

        simulated_data_file = open(dir, 'rw+')
        unique_rows = {}
        full_data_dict = {}
        full_info = []
        for line in simulated_data_file:
            line = line.split(" | ")
            cluster_label = int(line[1])
            vector = line[0][1:-1]
            full_info.append(self.DataNode(vector, cluster_label))

            if full_data_dict.get(cluster_label):
                full_data_dict[cluster_label] += [map(int, vector.split(","))]
            else:
                full_data_dict[cluster_label] = [map(int, vector.split(","))]
            if vector in unique_rows.keys():
                unique_rows[vector] += 1
            else:
                unique_rows[vector] = 1

        return unique_rows, full_data_dict, full_info

    '''
        Given single cell sequencing data, the method will process the file and would return similar output
        as the above method.
    '''

    def process_single_cell_data_file(self, dir):

        sample_data_file = open(dir, 'rw+')
        unique_rows = {}
        full_data_dict = {}
        full_info = []
        labels = (sample_data_file.readline()).split(",")
        data_lines = [[] for _ in range(len(labels) - 1)]
        for line in sample_data_file:

            line = (line.split('"')[2][1:])
            line = map(str, line.split(","))
            for i in xrange(len(line)):
                data_lines[i].append(line[i])

        for label, column in zip(labels[1:], data_lines):
            column = [x.replace("\n", "") for x in column]
            full_data_dict[label] = column
            row = ','.join(str(x) for x in column)
            if label != " ":
                full_info.append(self.DataNode(row, label))

            unique_rows[row] = 1

        return unique_rows, full_data_dict, full_info

    '''
        Get's averaged genotypes of clustered SCS data.
    '''

    def get_single_cell_genotypes(self, clustered_data_dict, full_or_data, delimiter=","):
        genotypes = []
        labels_in_clusters = {}
        for key, values in clustered_data_dict.items():
            genotype = [{'0': 0, '1': 0, 'NA': 0} for _ in range(len(values[0].split(delimiter)))]
            print "Cluster number:", key, " | Number of members:", len(values)
            for sample in values:
                for label, vector in full_or_data.items():
                    if vector == sample.split(delimiter):
                        if labels_in_clusters.get(key):
                            labels_in_clusters[key] += [label]
                        else:
                            labels_in_clusters[key] = [label]

                for index, el in enumerate(sample.split(delimiter)):
                    genotype[index][el] += 1
            genotype_list = []
            for index, el_dict in enumerate(genotype):
                genotype_list.append(max(el_dict, key=el_dict.get))
            genotypes.append(genotype_list)

        string_genotype = []
        for el in genotypes:
            string_genotype += [",".join(str(e) for e in el)]

        return string_genotype, labels_in_clusters

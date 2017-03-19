class Helper:
    def __init__(self):
        pass

    def get_label_of_cluster(self, vector, full_dict):
        vector = map(int, vector.split(","))
        for key in full_dict.keys():
            if vector in full_dict[key]:
                return key

    def diff_char_string(self, a, b):
        distance = 0
        for i, j in zip(a, b):
            if i != j:
                distance += 1
        return distance

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


    def read_simulated_data_file(self, dir):
        class DataNode:
            def __init__(self, vector, cluster_label):
                self.vector = vector
                self.cluster_label = cluster_label

        simulated_data_file = open(dir, 'rw+')
        unique_rows = {}
        full_data_dict = {}
        full_info = []
        for line in simulated_data_file:
            # line = line[:-1]
            # label = int(line.split(",")[-1])
            line = line.split(" | ")
            cluster_label = int(line[1])
            vector = line[0][1:-1]
            full_info.append(DataNode(vector, cluster_label))

            if full_data_dict.get(cluster_label):
                full_data_dict[cluster_label] += [map(int, vector.split(","))]
            else:
                full_data_dict[cluster_label] = [map(int, vector.split(","))]
            if vector in unique_rows.keys():
                unique_rows[vector] += 1
            else:
                unique_rows[vector] = 1

        return unique_rows, full_data_dict, full_info
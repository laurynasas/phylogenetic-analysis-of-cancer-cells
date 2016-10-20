import numpy as np


def process_file(dir):
	sample_data_file = open(dir, 'rw+')
	unique_rows = {}
	sample_data_file.readline()

	for line in sample_data_file:
		line = line[:-1]
		line = ','.join([str(x) for x in (line.split(",")[:-1])])

		if line in unique_rows.keys():
			unique_rows[line] +=1
		else:
			unique_rows[line] = 1

	return unique_rows


def diff_char_string(a,b):
	distance = 0
	for i,j in zip(a,b):
		if i != j:
			distance+=1
	if distance == 0:
		print a, b
	return distance


def find_distance_matrix(unique_rows):
	number_rows = len(unique_rows.keys())
	unique_strings = unique_rows.keys()

	distance_matrix = []
	for i in range(number_rows):
		new_row = []
		for j in range(number_rows):
			new_row +=[diff_char_string(unique_strings[i], unique_strings[j])]
		distance_matrix+=[new_row]

	return distance_matrix


def find_min_corner(distance_matrix):
	print distance_matrix
	minimal = max(distance_matrix[0][1:])
	location = (-1,-1)
	for i in range(len(distance_matrix[0])):
		for j in range(i,len(distance_matrix[0])):
			if distance_matrix[i][j] != 0 and distance_matrix[i][j] < minimal:
				minimal = distance_matrix[i][j]
				location = (i,j)
	print minimal, location
	return location


def delete_row_and_column(distance_matrix,location):
	for i in range(len(distance_matrix[0])):
		del distance_matrix[i][location[0]]
		del distance_matrix[i][location[1]-1]
	del distance_matrix[location[0]]
	del distance_matrix[location[1]-1]

	return distance_matrix


def add_new_row_and_column(distance_matrix,row, position):
	distance_matrix = distance_matrix[:position[0]-1] + [row] + distance_matrix[position[0]:]
	distance_matrix = distance_matrix[:position[1]-1] + [row] + distance_matrix[position[1]:]
	print np.matrix(distance_matrix)

def merge(distance_matrix, location, clusters):
	# print len(distance_matrix[0])	
	distance_matrix = distance_matrix.tolist()
	new_row = [ min(distance_matrix[location[0]][i],distance_matrix[location[1]][i]) for i in range(len(distance_matrix[0]))]
	del new_row[location[1]]
	print "new row: ",new_row
	distance_matrix = np.matrix(distance_matrix)

	clusters[location[0]] += [clusters[location[1]]]
	del clusters[location[1]]
	print clusters
	distance_matrix = np.delete(distance_matrix, (location[0]), axis=0)
	distance_matrix = np.delete(distance_matrix, (location[1]-1), axis=0)

	distance_matrix =np.delete(distance_matrix, (location[0]), axis=1)
	distance_matrix =np.delete(distance_matrix, (location[1]-1), axis=1)
	distance_matrix = np.asarray(distance_matrix)
	print "Before adding new lines\n",distance_matrix
	distance_matrix = np.insert(distance_matrix, location[0], new_row[:location[0]]+new_row[location[0]+1:], axis=0)
	print "Before adding new lines\n",distance_matrix

	distance_matrix = np.insert(distance_matrix, location[0], new_row, axis=1)
	distance_matrix = np.matrix(distance_matrix)
	print "After added new liens\n",distance_matrix

	return distance_matrix, clusters
	# print new_row
	# print distance_matrix
	print distance_matrix.size
	# delete_row_and_column(distance_matrix, location)
	# add_new_row_and_column(distance_matrix, new_row, location)

def single_linkage_clustering(distance_matrix):
	clusters = {}
	for x in range(distance_matrix[0].size):
		clusters[x] = [x]
	print clusters
	location = find_min_corner(distance_matrix.tolist())
	while location != (-1,-1):
		# print 
		location = find_min_corner(distance_matrix.tolist())
		if location != (-1,-1):
			(distance_matrix,clusters) = merge(distance_matrix, location, clusters)
	print "final",clusters
if __name__ == '__main__':
	dir = "./data/sample_data.csv"
	unique_rows = process_file(dir)
	print unique_rows
	# distance_matrix = np.matrix(find_distance_matrix(unique_rows))
	distance_matrix = np.matrix([[0,662,877,255,412,996],[662,0,295,468,268,400],[877,295,0,754,564,138],[255,468,754,0,219,869],[412,268,564,219,0,669],[996,400,138,869,669,0]])
	print "--->\n",distance_matrix
	print single_linkage_clustering(distance_matrix)
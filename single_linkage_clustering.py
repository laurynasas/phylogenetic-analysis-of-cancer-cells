import numpy as np


def find_min_corner(distance_matrix):
	# print distance_matrix
	minimal = max(distance_matrix[0][1:])
	location = (-1,-1)
	for i in range(len(distance_matrix[0])):
		for j in range(i,len(distance_matrix[0])):
			if distance_matrix[i][j] != 0 and distance_matrix[i][j] < minimal:
				minimal = distance_matrix[i][j]
				location = (i,j)
	print minimal, location
	return location


def merge(distance_matrix, location, clusters):
	# print len(distance_matrix[0])	
	distance_matrix = distance_matrix.tolist()
	new_row = [ min(distance_matrix[location[0]][i],distance_matrix[location[1]][i]) for i in range(len(distance_matrix[0]))]
	del new_row[location[1]]
	print "new row: ",new_row
	distance_matrix = np.matrix(distance_matrix)
	print "Cluster before failing", clusters
	clusters[location[0]] += [clusters[location[1]]]
	del clusters[location[1]]
	print "Cluster after adding new stuff",clusters
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
	print "SIZE after merging", distance_matrix.size

	return distance_matrix, clusters
	# print new_row
	# print distance_matrix
	# delete_row_and_column(distance_matrix, location)
	# add_new_row_and_column(distance_matrix, new_row, location)

def single_linkage_clustering(distance_matrix):
	clusters = []
	for i in range(distance_matrix[0].size):
		clusters += [['V'+str(i)]]
	print clusters
	location = find_min_corner(distance_matrix.tolist())
	while location != (-1,-1):
		location = find_min_corner(distance_matrix.tolist())
		if location != (-1,-1):
			(distance_matrix,clusters) = merge(distance_matrix, location, clusters)
	# print "final",clusters
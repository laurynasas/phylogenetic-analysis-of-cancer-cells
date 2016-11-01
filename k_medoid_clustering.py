import numpy as np
import random
from processing_sample_data import *
from sklearn.metrics.pairwise import pairwise_distances
from silhouette_score_implementation import *


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
		J = np.argmin(D[:,M], axis=1)
		for kappa in range(k):
			C[kappa] = np.where(J==kappa)[0]
		# update cluster medoids
		for kappa in range(k):
			# print D[np.ix_(C[kappa],C[kappa])]

			J = np.mean(D[np.ix_(C[kappa],C[kappa])],axis=1)
			if len(J)!=0:
				j = np.argmin(J)
				Mnew[kappa] = C[kappa][j]
		np.sort(Mnew)
		# check for convergence
		if np.array_equal(M, Mnew):
			break
		M = np.copy(Mnew)
	else:
		# final update of cluster memberships
		J = np.argmin(D[:,M], axis=1)
		for kappa in range(k):
			C[kappa] = np.where(J==kappa)[0]

	# return results
	return M, C

if __name__ == '__main__':
	dir = "./data/sample_data.csv"
	unique_rows = process_file(dir)
	distance_matrix = np.matrix(find_distance_matrix(unique_rows))
	# print unique_rows
	M,C = do_k_medoid_clustering(distance_matrix, 10)
	# print type(distance_matrix.shape[0])	
	labels = np.zeros(distance_matrix.shape[0])
	print C
	for key in zip(C.keys()):
		value_list = C[key[0]]
		for value in value_list:
			labels[value] = key[0]

	print silhouette_score_slow(distance_matrix,labels)
class Pipeline:
    def __init__(self, raw_data_dir, true_genotypes_dir, clustering_method, tree_method, number_of_clusters, vector_size):
        self.clustering_method = clustering_method
        self.raw_data_dir = raw_data_dir
        self.true_genotypes_dir = true_genotypes_dir
        self.no_clusters = number_of_clusters
        self.vector_size = vector_size
        self.tree_method = tree_method

    def run_pipe(self):

        if self.clustering_method == "slc":
            self._clustered_data_dict = self._do_slc()
        elif self.clustering_method == "k_means":
            self._clustered_data_dict = self._do_k_means()
        elif self.clustering_method == "bmm":
            self._clustered_data_dict = self._do_bmm()

        self.predicted_genotypes = self._get_non_unique_genotypes()
        self.true_genotypes = self._load_true_genotypes()

        self._equate_dimensions()

        self._save_genotypes()


        if self.tree_method == "slc":
            self._clustered_data_dict = self._do_slc()
        elif self.tree_method == "k_means":
            self._clustered_data_dict = self._do_k_means()
        elif self.tree_method == "bmm":
            self._clustered_data_dict = self._do_bmm()



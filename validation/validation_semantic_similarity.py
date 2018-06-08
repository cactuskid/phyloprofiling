import numpy as np

from utils import goatools_utils


class Validation_semantic_similarity(object):

    def __init__(self, go_file, term_counts, go_terms_hdf5):
        self.go_file = go_file
        self.term_counts = term_counts

        self.goterms2parents = go_terms_hdf5['goterms2parents']
        self.hog2goterms = go_terms_hdf5['hog2goterms']

    def semantic_similarity_score(self, hog_id_1, hog_id_2):
        """
        Runs semantic similarity analysis from 2 hog ids
        :param hog_id_1: first hog id
        :param hog_id_2: second hog id
        :return: semantic similarity score between the two hog ids
        """
        go_terms_1 = goatools_utils.get_go_terms_hdf5(hog_id_1, self.hog2goterms)
        go_terms_2 = goatools_utils.get_go_terms_hdf5(hog_id_2, self.hog2goterms)

        # if go_terms_1 and go_terms_2:
        score = self._compute_score(go_terms_1, go_terms_2)
        # else:
        #     score = -1
        return score

    def _compute_score(self, query_go_terms, result_go_terms):
        """
        Computes semantic similarity score between two hogs
        :param query_go_terms: dict of genes: list of go terms
        :param result_go_terms: dict of genes: list of go terms
        :return: semantic similarity score
        """
        # for each couple, compute resnik

        dist_mat = self._compute_genes_distance(query_go_terms, result_go_terms)

        # bma on this matrix

        score = self._mean_max_score_matrix(dist_mat)

        return score

    def _compute_genes_distance(self, go_terms_genes_1, go_terms_genes_2):
        """
        Computes matrix of distance between the genes of two hogs
        :param go_terms_genes_1: dictionary of genes
        :param go_terms_genes_2: dictionary of genes
        :return: matrix of distance between genes of hogs
        """
        # try:
        if type(go_terms_genes_1) is dict and type(go_terms_genes_2) is dict:
            keys_1 = go_terms_genes_1.keys()
            keys_2 = go_terms_genes_2.keys()

            gene_dist = np.zeros((len(keys_1), len(keys_2)))

            for k in range(len(keys_1)):
                for l in range(len(keys_2)):

                    go_terms_1 = list(go_terms_genes_1[list(keys_1)[k]])
                    go_terms_2 = list(go_terms_genes_2[list(keys_2)[l]])

                    if go_terms_1 and go_terms_2:
                        gene_dist[k, l] = self._compute_go_terms_score_per_gene(go_terms_1, go_terms_2)
            # except:
            #     gene_dist = -1
            return gene_dist
        else:
            return -1

    def _compute_go_terms_score_per_gene(self, go_terms_gene_1, go_terms_gene_2):
        """
        Computes the semantic similarity score between two genes
        :param go_terms_gene_1: list of go terms from one of the gene
        :param go_terms_gene_2: list of go terms from one of the gene
        :return: semantic similarity score between two genes
        """
        ss_dist = np.zeros((len(go_terms_gene_1), len(go_terms_gene_2)))

        for m in range(len(go_terms_gene_1)):
            for n in range(len(go_terms_gene_2)):
                dist = goatools_utils.resnik_sim_hdf5(go_terms_gene_1[m], go_terms_gene_2[n], self.go_file, self.term_counts, self.goterms2parents)
                ss_dist[m, n] = dist

        gene_score = self._mean_max_score_matrix(ss_dist)

        return gene_score

    @staticmethod
    def _mean_max_score_matrix(matrix):
        """
        Computes the BMA of a matrix
        :param matrix: matrix
        :return: score: BMA of the matrix; returns -1 if matrix has 0 or 1 dimension
        """
        #
        # try:
        matrix_size = np.prod(matrix.shape)
        if not matrix_size:
            return -1
        score = sum(matrix.max(0))+sum(matrix.max(1)) / matrix_size
        # except:
        #     score = -1

        return score

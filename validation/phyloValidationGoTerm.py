import numpy as np

from utils import goatools_utils

from time import time


class SemanticSimilarityAnalysis(object):

    def __init__(self, go, oma, term_counts, go_terms_hdf5):
        self.go_file = go
        self.oma = oma
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

        score = self._compute_score(go_terms_1, go_terms_2)

        return score

    def semantic_similarity_score_from_go_terms(self, go_terms_1, go_terms_2):
        """
        Runs semantic similarity analysis from 2 dictionaries of go terms
        :param go_terms_1: first dictionary of go terms
        :param go_terms_2: second dictionary of go terms
        :return: semantic similarity score between the two dictionaries of go terms
        """

        score = self._compute_score(go_terms_1, go_terms_2)

        return score

    def get_go_terms(self, hog_id):
        """
        Fetch the genes from hog id, then get all GO terms from it
        :param hog_id: hog id
        :return: list of GO term
        """
        return self._get_go_terms(hog_id)

    def _compute_score(self, query_go_terms, result_go_terms):
        """
        Computes semantic similarity score between two hogs
        :param query_go_terms: dict of genes: list of go terms
        :param result_go_terms: dict of genes: list of go terms
        :return: semantic similarity score
        """
        # for each couple, compute resnik
        start_time = time()

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
        keys_1 = go_terms_genes_1.keys()
        keys_2 = go_terms_genes_2.keys()

        gene_dist = np.zeros((len(keys_1), len(keys_2)))

        for k in range(len(keys_1)):
            for l in range(len(keys_2)):

                go_terms_1 = list(go_terms_genes_1[list(keys_1)[k]])
                go_terms_2 = list(go_terms_genes_2[list(keys_2)[l]])

                if go_terms_1 and go_terms_2:
                    gene_dist[k, l] = self._compute_go_terms_score_per_gene(go_terms_1, go_terms_2)

        return gene_dist

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
        matrix_size = np.prod(matrix.shape)
        if not matrix_size:
            return -1

        return sum(matrix.max(0))+sum(matrix.max(1)) / matrix_size

    def _get_go_terms(self, hog_id):
        """
        Fetch the genes from hog id, then get all GO terms from it
        :param hog_id: hog id
        :return: list of GO term
        """
        genes = self._get_hog_members(hog_id)
        go_dict = {entry: {self._format_go_term(e) for e in self._get_entry_gene_ontology(entry)} for entry in genes}
        go_dict_filtered = self._filter_result(go_dict)

        return self._clean_dictionary(go_dict_filtered)

    @staticmethod
    def _format_go_term(e):
        #return e['TermNr']
        return 'GO:{:07d}'.format(e['TermNr'])

    def _filter_result(self, go_dict):

        go_dict_filtered = {}

        for gene_name, terms in go_dict.items():
            filtered_terms = self._filter_namespace(terms)
            if filtered_terms:
                go_dict_filtered[gene_name] = filtered_terms

        return go_dict_filtered

    @staticmethod
    def _clean_dictionary(dictionary):
        return {k: v for k, v in dictionary.items() if v}

    def _get_entry_gene_ontology(self, entry):
        return self.oma.root.Annotations.GeneOntology.read_where('EntryNr == {}'.format(entry))

    def _get_hog_members(self, hog_id):
        """
        Gets all gene members from the hog
        :param hog_id: hog id
        :return: list of genes from the hog
        """
        iterator = self._iter_hog_member(hog_id)
        population = frozenset([x['EntryNr'] for x in iterator])
        return population

    @staticmethod
    def _hog_lex_range(hog):
        """
        Decodes hog format
        :param hog: (bytes or string): hog id
        :return: hog_str: encoded hog id
        """
        hog_str = hog.decode() if isinstance(hog, bytes) else hog
        return hog_str.encode('ascii'), (hog_str[0:-1] + chr(1 + ord(hog_str[-1]))).encode('ascii')

    def _iter_hog_member(self, hog_id):
        """
        iterator over hog members
        :param hog_id: hog id
        :return: yields members of hog
        """
        hog_range = self._hog_lex_range(hog_id)
        it = self.oma.root.Protein.Entries.where('({!r} <= OmaHOG) & (OmaHOG < {!r})'.format(*hog_range))
        for row in it:
            yield row.fetch_all_fields()

    def _filter_namespace(self, list_terms, name='biological_process'):
        """
        Keep only go terms within the correct ontology
        :param list_terms: list of go terms
        :param name: namespace to keep; default: 'biological_process'
        :return: list of terms with the correct namespace
        """
        terms_to_keep = []
        for term in list_terms:
            try:
                if self.go_file[term].namespace == name:
                    terms_to_keep.append(term)
            except KeyError:
                pass
        return terms_to_keep

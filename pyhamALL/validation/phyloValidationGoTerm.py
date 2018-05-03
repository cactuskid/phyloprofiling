
from goatools.semantic import resnik_sim
import numpy as np

"""
Created on Sun Mar 25 11:56:43 2018

@author: Laurent Kilchoer
"""
'''
This script takes a list of HOG ids and returns a matrix with all the semantic
distance between them.

The semantic distance score is computed using a BMA-RESNIK-BMA strategy
'''


class SemanticSimilarityAnalysis(object):

    def __init__(self, go, h5file, termcounts):
        self.go = go
        self.h5file = h5file
        self.termcounts = termcounts

    def semantic_similarity_score(self, query, result):
        ''' Runs semantic similarity analysis
        Args:
            query: query hog id
            result: result hog ids
        Returns:
            score : semantic similarity score between query hog and result hog
        '''
        query_go_terms = self.get_go_terms(query)
        result_go_terms = self.get_go_terms(result)

        score = self.compute_score(query_go_terms, result_go_terms)

        return score

    def compute_score(self, query_go_terms, result_go_terms):
        '''
        Computes semantic similarity score between two hogs
        Args:
            query_go_terms: dict of genes:list of go terms
            result_go_terms: dict of genes:list of go terms
        Returns:
            score: semantic similarity score
        '''
        # 1. compute distances between all genes
        distance_mat = self.compute_genes_distances(query_go_terms, result_go_terms)
        # 2. from all the distance, select BMA
        score = self.mean_max_score_matrix(distance_mat)

        return score

    def compute_genes_distances(self, query_go_terms, result_go_terms):
        '''
        Compute matrix of distances between the genes of two hogs
        Args:
            query_go_terms: dict of genes:list of go terms
            result_go_terms: dict of genes:list of go terms
        Returns:
            gene_dist: matrix of distance between genes of hogs
        '''
        query_keys = query_go_terms.keys()
        result_keys = result_go_terms.keys()

        gene_dist = np.zeros((len(query_keys), len(result_keys)))

        for k in range(len(query_keys)):
            for l in range(len(result_keys)):

                goterms1 = list(query_go_terms[list(query_keys)[k]])
                goterms2 = list(result_go_terms[list(result_keys)[l]])

                if goterms1 and goterms2:
                    gene_dist[k, l] = self.compute_genes_score(goterms1, goterms2)

        return gene_dist

    def compute_genes_score(self, goterms1, goterms2):
        '''
        Computes the semantic similarity score between two genes
        Args:
            goterms1: list of go terms from one of the gene of the query
            goterms2: list of go terms from one of the gene of the result
        Returns:
            gene_score: semantic similarity score between two genes
        '''
        ss_dist = np.zeros((len(goterms1), len(goterms2)))

        for m in range(len(goterms1)):
            for n in range(len(goterms2)):

                try:
                    dist = resnik_sim(goterms1[m], goterms2[n], self.go, self.termcounts)
                    ss_dist[m, n] = dist
                except:
                    # TODO catch real error
                    pass
        gene_score = self.mean_max_score_matrix(ss_dist)

        return gene_score

    def mean_max_score_matrix(self, sem_mat):
        '''
        computes the BMA of a matrix
        Args:
            sem_mat: scores matrix
        Returns:
            score: BMA of the scores matrix
        '''
        score = -1.
        for row in sem_mat:
            try:
                if score < 0:
                    score = 0
                score += max(i for i in row if i >= 0)
            except:
                # TODO catch actual error
                pass
        for col in sem_mat.transpose():
            try:
                if score < 0:
                    score = 0
                score += max(i for i in col if i >= 0)
            except:
                # TODO catch actual error
                pass
        mat_size = sum(sem_mat.shape)
        score /= mat_size

        return score

    # returns go terms
    def get_go_terms(self, hog_id):
        '''Fetch the genes from hog id, then get all GO terms from it
        Args:
            hog_id: hog id
        Returns:
            golist: list of GO term
        '''
        population = self.get_hog_members(hog_id)
        # turn into godict
        # for each gene, we have the go terms
        godict = {entry: {('GO:{:07d}'.format(e['TermNr'])) for e in self._get_entry_geneOntology(entry)} for entry in population}

        godict = {k: v for k, v in godict.items() if v}
        godictfinal = {}
        for gene, terms in godict.items():
            newTerms = self.filter_namespace(terms)
            if newTerms:
                godictfinal[gene] = newTerms
        godictfinal = {k: v for k, v in godictfinal.items() if v}
        return godictfinal

    def _get_entry_geneOntology(self, entry):
        return self.h5file.root.Annotations.GeneOntology.read_where('EntryNr == {}'.format(entry))

    def get_hog_members(self, hog_id, maxEntries=None):
        '''Get all gene members from the hog
        Args:
            hog_id: hog id
            maxEntries: max entries before breaking
        Returns:
            population: list containing the genes of the Hog
        '''
        iterator = self.iter_hog_member(hog_id)

        if maxEntries is None:
            population = frozenset([x['EntryNr'] for x in iterator])
        else:
            population = []
            for i, x in enumerate(iterator):
                if i > maxEntries:
                    break

                population.append(x['EntryNr'])
            population = list(population)

        return population

    # decode hog format
    def _hog_lex_range(self, hog):
        '''decode hog format
        Args:
            hog (bytes or string): hog id
        Returns:
            hog_str: encoded hog id
        '''
        hog_str = hog.decode() if isinstance(hog, bytes) else hog
        return hog_str.encode('ascii'), (hog_str[0:-1] + chr(1 + ord(hog_str[-1]))).encode('ascii')

    # yield hog members given hog id and oma database
    def iter_hog_member(self, hog_id):
        '''iter over hog members
        Args:
            hog_id: hog id
        Returns:
            yield members of hog
        '''
        hog_range = self._hog_lex_range(hog_id)
        it = self.h5file.root.Protein.Entries.where('({!r} <= OmaHOG) & (OmaHOG < {!r})'.format(*hog_range))
        for row in it:
            yield row.fetch_all_fields()

    def filter_namespace(self, list_terms, name='biological_process'):
        '''
        return list of terms with the correct namespace
        Args:
            list_terms: list of go terms
            name: namespace to keep
        Returns:
            terms_to_keep: list of terms with the correct namespace
        '''
        terms_to_keep = []
        for term in list_terms:
            try:
                if self.go[term].namespace == name:
                    terms_to_keep.append(term)
            except KeyError:
                pass
        return terms_to_keep

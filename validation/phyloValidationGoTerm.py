from goatools.semantic import resnik_sim
import numpy as np

from time import time

class SemanticSimilarityAnalysis(object):

    def __init__(self, go, h5file, termcounts):
        self.go = go
        self.h5file = h5file
        self.termcounts = termcounts
        self.used_hog_dict = {}

    def semantic_similarity_score(self, query, result):
        ''' Runs semantic similarity analysis
        Args:
            query: query hog id
            result: resulquery_go_termst hog ids
        Returns:
            score : semantic similarity score between query hog and result hog
        '''
        from time import time
        start_time = time()

        if query not in self.used_hog_dict:
            query_go_terms = self.get_deepest_go_term_per_gene(self.get_go_terms(query))
            self.used_hog_dict[query] = query_go_terms
        else:
            query_go_terms = self.used_hog_dict[query]

        if result not in self.used_hog_dict:
            result_go_terms = self.get_deepest_go_term_per_gene(self.get_go_terms(result))
            self.used_hog_dict[result] = result_go_terms
        else:
            result_go_terms = self.used_hog_dict[result]


        #score1 = self.compute_score(query_go_terms, result_go_terms)

        score = self.compute_score_method_2(query_go_terms, result_go_terms)
        #print("score 1 is {}, score 2 is {}".format(score1, score))
        return score

    def API_get_best_go_terms_per_gene(self, hog_id):
        go_terms = self.get_go_terms(hog_id)
        best_go_terms = self.get_deepest_go_term_per_gene(go_terms)
        return best_go_terms

    def API_compute_score_method_2(self, go_terms_1, go_terms_2):
        dist_mat = self.compute_genes_distances(go_terms_1, go_terms_2)
        score = self.mean_max_score_matrix(dist_mat)
        return score


    def compute_score_method_2(self, query_go_terms, result_go_terms):
        # new function !

        # select the best go per gene, and change the dict

        #deepest_query_go_terms = self.get_deepest_go_term_per_gene(query_go_terms)
        #deepest_result_go_terms = self.get_deepest_go_term_per_gene(result_go_terms)

        # for each couple, compute resnik
        dist_mat = self.compute_genes_distances_list_go_terms(query_go_terms, result_go_terms)
        # bma on this matrix
        score = self.mean_max_score_matrix(dist_mat)

        return score

    def get_deepest_go_term_per_gene(self, go_terms_dict):
        # new function !

        sign_go_terms = []
        for gen_name, go_terms in go_terms_dict.items():
            sign_go_terms.append(self.get_deepest_go_term(go_terms))

        sign_go_terms = list(set(sign_go_terms))
        return sign_go_terms

    def get_deepest_go_term(self, go_terms):
        return max(go_terms, key=lambda t: self.go[t].depth)

    def compute_genes_distances_list_go_terms(self, query_go_terms, result_go_terms):
        # new function !

        gene_dist = np.zeros((len(query_go_terms), len(result_go_terms)))

        for k in range(len(query_go_terms)):
            for l in range(len(result_go_terms)):
                try:
                    gene_dist[k, l] = resnik_sim(query_go_terms[k], result_go_terms[l], self.go, self.termcounts)
                except:
                    pass
        return gene_dist

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

        mat_size = np.prod(sem_mat.shape)
        if not mat_size:
            return -1

        return sum(sem_mat.max(0))+sum(sem_mat.max(1)) / mat_size

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
        start_time = time()
        godict = {entry: {('GO:{:07d}'.format(e['TermNr'])) for e in self._get_entry_geneOntology(entry)} for entry in population}
        print("get go terms {}".format(time()-start_time))
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
        start_time = time()
        if maxEntries is None:
            population = frozenset([x['EntryNr'] for x in iterator])
        else:
            population = []
            for i, x in enumerate(iterator):
                if i > maxEntries:
                    break

                population.append(x['EntryNr'])
            population = list(population)
        print('time for hog members stuff {}'.format(time()-start_time))
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

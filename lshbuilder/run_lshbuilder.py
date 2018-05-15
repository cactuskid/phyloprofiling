from tables import *
import functools
import multiprocessing as mp

from pyoma.browser import db

from utils import config_utils
from utils import files_utils
from utils import pyhamutils
from utils import hashutils

from lshbuilder import lsh_builder_utils


if __name__ == '__main__':

    with open_file(config_utils.omadir + 'OmaServer.h5', mode="r") as h5OMA:
        db_obj = db.Database(h5OMA)
        oma_id_obj = db.OmaIdMapper(db_obj)
        dic, tree = files_utils.get_species_tree_replacement_dic(h5OMA, oma_id_obj)
        taxaIndex, reverse = files_utils.generate_taxa_index(tree)

        HAM_PIPELINE = functools.partial(pyhamutils.get_ham_treemap, species_tree=tree)
        HASH_PIPELINE = functools.partial(hashutils.tree2hashes_from_row,
                                          events=['duplication', 'gain', 'loss', 'presence'], combination=True)
        ROW_PIPELINE = functools.partial(hashutils.tree2mat, taxaIndex=taxaIndex)

        columns = len(taxaIndex)
        rows = len(h5OMA.root.OrthoXML.Index)

        lsh_builder_utils.mp_with_timeout(number_workers=int(mp.cpu_count() / 2),
                                          number_updaters=1,
                                          data_generator=lsh_builder_utils.generates_dataframes(100),
                                          worker_function=lsh_builder_utils.worker,
                                          update_function=lsh_builder_utils.saver)

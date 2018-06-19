import profiler
from utils import config_utils
from datetime import datetime

date_string = "{:%B_%d_%Y_%H_%M}".format(datetime.now())
profilerObj = profiler.Profiler(lsh_path=config_utils.datadirLaurent + 'June_05_2018_14_31_0.7_newlsh.pkl',
                                hashes_path=config_utils.datadirLaurent + 'June_05_2018_14_31hashes.h5',
                                obo_file_path=config_utils.datadirLaurent + 'project/data/go.obo',
                                gaf_file_path=config_utils.datadirLaurent + 'project/data/gene_association.tair',
                                h5_go_terms_parents_path=config_utils.datadirLaurent + 'project/data/parents.h5')
date_string = "{:%B_%d_%Y_%H_%M}".format(datetime.now())
profilerObj.validate_pipeline(path_to_hog_id_file=config_utils.datadirLaurent + 'project/results/hog_list.csv',
                              path_to_save=config_utils.datadirLaurent + 'validation' + date_string + '.txt')

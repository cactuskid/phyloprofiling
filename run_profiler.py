import profiler
from utils import config_utils

if __name__ == '__main__':

    profiler = profiler.Profiler(config_utils.omadir + 'OmaServer.h5')
    profiler.go_benchmarking_init(config_utils.datadirLaurent + 'project/data/go.obo', config_utils.datadirLaurent + 'project/data/gene_association.tair', config_utils.datadirLaurent + 'project/data/parents.h5',)

    profiler.lsh_loader(config_utils.datadirLaurent + 'May_16_2018_17_33_0.7_newlsh.pkl')

    profiler.results_save(fam_id=789, scores=True, path_to_save=config_utils.datadirLaurent + 'test_01')


%load_ext autoreload
%autoreload 2
import profiler
from utils import config_utils
profilerObj = profiler.Profiler(config_utils.omadir + 'OmaServer.h5')
profilerObj.go_benchmarking_init(config_utils.datadirLaurent + 'project/data/go.obo', config_utils.datadirLaurent + 'project/data/gene_association.tair', config_utils.datadirLaurent + 'project/data/parents.h5')
profilerObj.lsh_loader(config_utils.datadirLaurent + 'May_23_2018_18_29_0.7_newlsh.pkl', config_utils.datadirLaurent + 'May_23_2018_18_29hashes.h5')


profilerObj.save_results(fam_id=589123, path_to_save=config_utils.datadirLaurent + 'test_06.txt')



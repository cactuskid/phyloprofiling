import profiler
from utils import config_utils

if __name__ == '__main__':

    profiler = profiler.Profiler(config_utils.omadir + 'OmaServer.h5')
    profiler.go_benchmarking_init(config_utils.datadirLaurent + 'project/data/go.obo', config_utils.datadirLaurent + 'project/data/gene_association.tair')

    profiler.lsh_loader(config_utils.datadirLaurent + 'May_16_2018_17_33_0.7_newlsh.pkl')

    profiler.results_save(fam_id=789, scores=True, path_to_save=config_utils.datadirLaurent + 'test_01')


%load_ext autoreload
%autoreload 2
import profiler
from utils import config_utils
profiler = profiler.Profiler(config_utils.omadir + 'OmaServer.h5')
profiler.go_benchmarking_init(config_utils.datadirLaurent + 'project/data/go.obo', config_utils.datadirLaurent + 'project/data/gene_association.tair')
profiler.lsh_loader(config_utils.datadirLaurent + 'May_16_2018_17_33_0.7_newlsh.pkl')
profiler.results_save(fam_id=789, scores=True, path_to_save=config_utils.datadirLaurent + 'test_01')
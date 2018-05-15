import profiler

if __name__ == '__main__':

    profiler = profiler.Profiler('path to oma')
    profiler.go_benchmarking_init('obo file path', 'gaf file path')

    profiler.lsh_loader('lsh file path')

    profiler.results_save(hog_id= 789, scores=True, path_to_save='path to save')

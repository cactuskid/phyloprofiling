
#!/bin/bash
#BSUB -L /bin/bash

#BSUB -o BOout2.txt -o:	standard	output
#BSUB -e BOerror2.txt

pyenv activate myenv
python lshoptimizer.py  --dir /scratch/cluster/monthly/dmoi/profiling/run_all2/ --db all

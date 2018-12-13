
#!/bin/bash
#BSUB -L /bin/bash

#BSUB -o BOout2.txt
#BSUB -e BOerror2.txt
#BSUB -M 100000000
#BSUB -R "rusage[mem=100000]"
#BSUB  -W 1000

source /home/dmoi/.pyenv/versions/myenv/bin/activate
python lshoptimizer.py  --pw 1  --lw 0  --dw 0 --dir /scratch/cluster/monthly/dmoi/profiling/run_all2/ --db all

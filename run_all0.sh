
#!/bin/bash
#BSUB -L /bin/bash

#BSUB -o BOout0.txt
#BSUB -e BOerror0.txt
#BSUB -M 100000000
#BSUB -R "rusage[mem=100000]"
#BSUB  -W 1000

source /home/dmoi/.pyenv/versions/myenv/bin/activate
python lshoptimizer.py  --pw 0  --lw 0  --dw 1 --dir /scratch/cluster/monthly/dmoi/profiling/run_all0/ --db all

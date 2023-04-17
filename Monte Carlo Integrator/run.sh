#!/bin/bash

#SBATCH ...

export OMP_STACKSIZE=16000M
export OMP_NUM_THREADS=64

module load buildenv/default-intel-2022a
module load SciPy-bundle/2022.05-intel-2022a 


python3

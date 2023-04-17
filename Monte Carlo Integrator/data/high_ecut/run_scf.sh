#!/bin/bash

#SBATCH ...

module load buildenv/default-intel-2022a
export KMP_NUM_THREADS=32
export KMP_STACKSIZE=7500M

../../../q-e-qe-6.4.1/PW/src/pw.x <scf.in> scf.out


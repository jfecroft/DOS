#!/bin/bash
#PBS -V
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=048:00:00
#PBS -o /$PBS_O_WORKDIR/$PBS_JOBID.stdout
#PBS -e /$PBS_O_WORKDIR/$PBS_JOBID.stderr
#PBS -N vibrations	


cd $PBS_O_WORKDIR
python rovib.py

#!/bin/bash -l
#PBS -A YOUR_PROJECT_ID
#PBS -N postproc_flat
#PBS -o logs/postproc_flat.log
#PBS -e logs/postproc_flat.log
#PBS -l walltime=24:00:00
#PBS -q YOUR_QUEUE
#PBS -l select=1:ncpus=18:mem=1400GB:ngpus=0
## preempt=0.2, economy=0.7, regular=1, premium=1.5
##PBS -l job_priority=premium
#PBS -M YOUR_EMAIL
#PBS -m abe
#PBS -r n

# Clear the environment from any previously loaded modules
module purge
module load ncarenv/25.10 gcc ncarcompilers netcdf
module li

#/glade/u/apps/ch/opt/usr/bin/dumpenv # Dumps environment (for debugging with CISL support)

time python 00_postproc_flat.py 2>&1 | tee logs/postproc_flat.out

qstat -f $PBS_JOBID >> logs/postproc_flat.out

#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --job-name={simname_ascii}
#SBATCH --output=logs/{simname_ascii}.log
#SBATCH --error=logs/{simname_ascii}.log
#SBATCH --qos=regular
#SBATCH --gpus-per-node=1
#SBATCH --mem=48GB
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=YOUR_EMAIL
{options_string}

export SLURM_CPU_BIND="cores"

echo $CUDA_VISIBLE_DEVICES

cd $SCRATCH/turbulent-seamount/simulations/  # or your clone path
juliaup default 1.12
time julia --project --pkgimages=no {julia_script} {run_options} --simname={simname} 2>&1 | tee logs/{simname_ascii}.out

#scontrol show job $PBS_JOBID >> logs/{simname_ascii}.log
#scontrol show job $PBS_JOBID >> logs/{simname_ascii}.out

#!/bin/bash
#SBATCH --job-name=prova_parallel
#SBATCH --output logs/%x.%j.out
#SBATCH --time=00:30:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=1

echo "SLURM_NODELIST=$SLURM_NODELIST"
echo "SLURM_NTASKS=$SLURM_NTASKS"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"


python3 simulation_bora.py

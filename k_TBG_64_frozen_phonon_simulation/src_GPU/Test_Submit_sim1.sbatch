#!/bin/bash
#SBATCH --job-name=s1-1
#SBATCH -A m5007_g
#SBATCH -C gpu
#SBATCH -t 0:29:50
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --gpus-per-task=1
#SBATCH --qos debug

# module load python
 conda activate abTEM

mexcuda -lcufft -R2018a cal_grad_all.cu
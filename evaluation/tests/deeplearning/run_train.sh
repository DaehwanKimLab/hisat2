#!/bin/bash -l
#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=GPUv1
#SBATCH --gres=gpu:1
#SBATCH --workdir=/home2/s175443/work/image_analysis/leela-zero/training/tf
/home2/s175443/work/image_analysis/leela-zero/training/tf/hisatdl_train.py

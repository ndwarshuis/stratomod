#!/bin/bash
#SBATCH --job-name=train_indel_ebm
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --mem=64G

cd /aigenomics/EBM_dev
module load anaconda
conda activate interpretML_env

python running_indels_ebm_training.py

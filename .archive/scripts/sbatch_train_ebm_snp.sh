#!/bin/bash
#SBATCH --job-name=train_snp_ebm
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --mem=64G

cd /aigenomics/EBM_dev
module load anaconda
conda activate interpretML_env

python running_snps_ebm_training.py

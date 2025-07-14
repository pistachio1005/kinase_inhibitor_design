#!/usr/bin/env bash
#SBATCH --job-name=pkn2_caga_kstar
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160gb
#SBATCH --time=96:00:00
#SBATCH --partition=compsci  
#SBATCH --output=/home/users/ys472/ying_Project/kinase_inhibitor_design/exelogs/%A_%a.out
#SBATCH --error=/home/users/ys472/ying_Project/kinase_inhibitor_design/exelogs/%A_%a.err

# Activate the conda base environment
source /home/users/ys472/miniconda3/bin/activate AmberTools22

# Run the Python script
python3 /home/users/ys472/ying_Project/kinase_inhibitor_design/MFS_run_test.py

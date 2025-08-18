#!/bin/csh -f
#SBATCH --output=submit.out
#SBATCH --mem=750G
#SBATCH --cpus-per-task=48
#SBATCH -A grisman --partition=grisman
python3 ccKStar_csv.py ./complex.ccsx ./target.ccsx ./design.ccsx
rm *confdb
echo "completed"
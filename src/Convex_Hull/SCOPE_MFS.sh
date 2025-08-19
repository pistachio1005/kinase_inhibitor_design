#!/bin/csh -f
#SBATCH --output=mfs.submit.out
#SBATCH --error=mfs.submit.err
#SBATCH --mem=256GB
#SBATCH --cpus-per-task=4
#SBATCH --partition=biostat

# Note: Dropped #SBATCH --cpus-per-task=48, changed to 4

set FILENAME = $1
set ARG = "$FILENAME"
set CMD = "python3 SCOPE_MFS.py $ARG"
eval $CMD
echo "Completed!"
#!/bin/csh -f
#SBATCH --output=submit.out
#SBATCH --error=submit.err
#SBATCH --mem=750GB
#SBATCH --cpus-per-task=48
#SBATCH --partition=biostat

#set FILENAME = $1
#set ARG = "$FILENAME"
set CMD = "python3 MFS_run.py"# $ARG"

eval $CMD

echo "Completed!"

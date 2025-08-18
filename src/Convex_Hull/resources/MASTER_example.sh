#!/bin/csh -f
#SBATCH --output=submit.out
#SBATCH --mem=750G
#SBATCH --cpus-per-task=48
#SBATCH --partition=grisman
./master --query 8gal-renamed-inverted-design.pds --targetList ~/dlab/henry/master-db/db.txt --rmsdCut 10.0 --topN 40 --outType match --seqOut ./matches-8gal.txt --structOut ./matches-8gal
echo "completed"
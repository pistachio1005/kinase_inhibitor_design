# example main file for D-space Thanatin designs. This file can be repurposed for diverse design objectives.
import glob
import os
import shutil

# -------------------------- MASTER --------------------------

"""
Recommended workflow for D-space designs:
1. rename the protonated L-L pdb using complex_renamer (target to A, design to z)
2. create MASTER D-space query- invert and isolate the L-design using reflect_design
3. STOP: run MASTER (code not in this project) on the D-design to get matches
4. rename all MASTER returns to chain B using chain_relabel
5. protonate the L-space MASTER matches (AMBERTOOLS can't protonate D-space) - we need bb protons for hull alignment
6. generate D-peptide:L-protein compelxes using scaffold_generator
"""

from MASTERPrep import *

# specify the L:L pdb (protonated) to rename
# recommended to use A/z naming to avoid issues with PDBs having identical chain names
# this code assumes the target ID < design ID (e.g., target B and design A will fail)
# complex_renamer('8gal.pdb', 'A', 'z')

# get the D-space design chain using new chain names
# reflect_design('8gal-renamed.pdb', 'A', 'z')

# STOP: now run MASTER (code not in this repo); see resources/MASTER_example.sh for an example script

# rename MASTER returns (in-place) to B (avoid name clash)
# chain_relabel('matches-8gal', 'B')

# add protons to these L-space backbones and the L-target (for later use for CH placement)
# We haven't corrected atom labels so some sc protons will be wrong, but this is ok. We only care about bb protons.
# protonate_chains('matches-8gal', '8gal-renamed.pdb')

# reflect MASTER returns to D-space and generate D-design:L-protein PDBs
# scaffold_generator("8gal-renamed.pdb", 'z', "matches-8gal", "scaffolds")

# STOP: I recommend reviewing the complexes in scaffolds to ensure correctness
# check that design residues are D, target is L, correct protonation on backbone for both chains, correct docking
# severe clashes are to be expected due to MASTER alignment. We'll fix this with ABAA.

# -------------------------- ABAA --------------------------

"""
Affinity-based Backbone Alignment Adjustment (ABAA)
The purpose of this algorithm is to update scaffold components (D-peptide and L-target) to their new chemical context by
using an affinity-based search (K*) to translate/rotate the design backbone while flexing target residues

Algorithm Overview (for each scaffold):
1. All D-design residues are reassigned to ALA
2. A VAL hull (approximately an ALA residue + VdW radius) is placed on each design residue
3. A flexibility hull (hull of WT identity rotamers) is placed on each L-target residue
4. Intra-chain hull intersection = flexibility on participating target residue
5. K* algorithm is run with flexible target residue and translation+rotation of design chain
6. GMEC structure of K* ensemble output is now adjusted and viable for Pair Set K* Search
"""

from KStarPrep import *

# for file in os.listdir("scaffolds"):
#
#     print("\n\n--- Now running ABAA on %s ---\n\n" % file)
#
#     # set up the output directory
#     full_filepath = os.path.join("scaffolds", file)
#     shutil.copy(full_filepath, file)
#     new_foldername = ("%s-ABAA" % (file.split("-")[0]))
#     os.mkdir(new_foldername)
#
#     # find contacts and prepare OSPREY K* files
#     # maintains GLY and PRO residues, which are likely important for bb conformation
#     osprey_fileprep_ABAA(file, new_foldername, 'B', 'D', False, True, "resources/K_bash.sh")
#     os.remove(file)

# # organize into a single directory
# print("\n\nMoving all matches into ABAA directory\n\n")
# os.mkdir("ABAA")
# for file in glob.glob("*-ABAA"):
#     final_path = os.path.join("ABAA", file)
#     shutil.move(file, final_path)

# STOP: now run these file on the cluster. Use resources/abaa_cluster_runner.sh for fast queueing
# recommendation: use resources/get_ABAA_GMEC.py for fast GMEC file organization (after runs are complete)


# -------------------------- Full Sequence Verification --------------------------

"""
Full Sequence Verification
The purpose of Full Sequence Verification (FSV) is to dynamically assign residues mutations to construct
full sequences.
"""

# rec: use FSV.sh to run FSV on the cluster



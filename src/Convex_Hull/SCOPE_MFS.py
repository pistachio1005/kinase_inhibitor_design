import sys

from pdbfixer import PDBFixer
from Find_Doublets import singlechain_design_info

import osprey
osprey.start()

# Parse command line arguments
if len(sys.argv) < 2:
    print("Usage: python SCOPE_MFS.py <file_name>")
    sys.exit(1)

file_name = sys.argv[1]

# Search for possible clashes with SCOPE
fixer = PDBFixer(filename=file_name)

residue_count_chain_b = 0
for chain in fixer.topology.chains(): #
    if chain.id == 'B':
        for residue in chain.residues(): #
            residue_count_chain_b += 1

print(f"Number of residues in Chain B: {residue_count_chain_b}")

# Parameterization
pdb_file = file_name
output_folder = "output_dir" # Output directory for hull PDBs
design_chain = "B"                    
designable_aas = ["TRP"]# ["VAL", "LEU"]       
save_pdbs = True # Check this...
chirality = "L" # "L"|"D"
fixed_residues = [i for i in range (1, residue_count_chain_b+1)] # This lets us keep residues WT (I think)

print("Beginning SCOPE.")

# Define and Run
contacts, interchain = singlechain_design_info(
    pdb_file,
    output_folder,
    design_chain,
    designable_aas,
    save_pdbs,
    chirality,
    fixed_residues
)

print("Completed SCOPE.")

# Pull important values

# Internal (ligand) clashes
int_clashes = [x for s in contacts if len(s) > 1 for x in s]

# inter-chain clashes
flat = [x for sub in interchain for x in sub]
seen = set()
duplicates = set()

for num in flat:
    if num in seen:
        duplicates.add(num)
    else:
        seen.add(num)
        
# Pull indices with more than average clashes
avg_len = sum(len(sub) for sub in interchain) / len(interchain)   # average number of ints per sublist
indices = [i for i, sub in enumerate(interchain) if len(sub) >= avg_len]
ligand_residues = [i+1 for i in indices].append(int_clashes)

# Flatten the list and remove duplicates using a set
flat = []
if ligand_residues:
    for item in ligand_residues:
        if isinstance(item, list):
            flat.extend(item)
        else:
            flat.append(item)

# Remove duplicates while preserving order
ligand_result = list(dict.fromkeys(flat))

# Protein
print("Making the following residues on protein flexible:")
print(duplicates)
print("Making the following residues on ligand flexible:")
print(ligand_result)

print("Completed SCOPE!")

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# read a PDB file for molecular info
mol = osprey.readPdb(file_name)

# define the protein strand-- here we want the entire protein(?)
protein = osprey.Strand(mol, templateLib=templateLib, residues = ['A1', 'A311'])

flex_residues = ["A"+str(i) for i in duplicates]

for res in flex_residues:
    protein.flexibility[res].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# Define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['B1', 'B10'])#, residues=['E948', 'E961'])

# These residues (3, 5, 9, 12), if changed, result in decreased binding capacity of MKI
flex_residues = ["B"+str(i) for i in ligand_result]

for res in flex_residues:
    # Drop the wild-type residue from the list we're searching over
    ligand.flexibility[res].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# Set translation/rotation
trflex = osprey.c.confspace.StrandFlex.TranslateRotate(10, 2.0)

ligandForConf = [ligand, trflex]

# make the conf space for the protein
proteinConfSpace = osprey.ConfSpace([protein])

# make the conf space for the ligand
ligandConfSpace = osprey.ConfSpace([ligandForConf])

# make the conf space for the protein+ligand complex
complexConfSpace = osprey.ConfSpace([protein, ligandForConf])

# how should we compute energies of molecules?
# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism)

# configure K*
kstar = osprey.KStar(
	proteinConfSpace,
	ligandConfSpace,
	complexConfSpace,
	epsilon=0.683, # you proabably want something more precise in your real designs
	writeSequencesToFile='MFS.results19Aug25.tsv'#,
    #maxSimultaneousMutations=100 # Added this to see if we could get the entire molecule converted to alanine
)

# configure K* inputs for each conf space
for info in kstar.confSpaceInfos():
	# how should we define energies of conformations?
	eref = osprey.ReferenceEnergies(info.confSpace, ecalc)
	info.confEcalc = osprey.ConfEnergyCalculator(info.confSpace, ecalc, referenceEnergies=eref)

	# compute the energy matrix
	emat = osprey.EnergyMatrix(info.confEcalc, cacheFile='emat.%s.dat' % info.id)

	# how should we score each sequence?
	# (since we're in a loop, need capture variables above by using defaulted arguments)
	def makePfunc(rcs, confEcalc=info.confEcalc, emat=emat):
		return osprey.PartitionFunction(
			confEcalc,
			osprey.AStarTraditional(emat, rcs, showProgress=False),
			osprey.AStarTraditional(emat, rcs, showProgress=False),
			rcs
		)
	info.pfuncFactory = osprey.KStar.PfuncFactory(makePfunc)

# run K*
scoredSequences = kstar.run(ecalc.tasks)

# make a sequence analyzer to look at the results
analyzer = osprey.SequenceAnalyzer(kstar)

# use results
for scoredSequence in scoredSequences:
	print("result:")
	print("\tsequence: %s" % scoredSequence.sequence)
	print("\tK* score: %s" % scoredSequence.score)

	# write the sequence ensemble, with up to 10 of the lowest-energy conformations
	numConfs = 10
	analysis = analyzer.analyze(scoredSequence.sequence, numConfs)
	print(analysis)
	analysis.writePdb(
		'seq.%s.pdb' % scoredSequence.sequence,
		'Top %d conformations for sequence %s' % (numConfs, scoredSequence.sequence)
	)
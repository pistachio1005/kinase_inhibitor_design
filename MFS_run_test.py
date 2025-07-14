import osprey
import os
osprey.start(heapSizeMiB=24000)

# Set up output directory
output_dir = "kstar_results"
os.makedirs(output_dir, exist_ok=True)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# read a PDB file for molecular info
mol = osprey.readPdb('/home/users/ys472/ying_Project/kinase_inhibitor_design/pkn2.cleaned.filtered.pdb')

# define the protein strand-- here we want the entire protein(?)
protein = osprey.Strand(mol, templateLib=templateLib, residues=['A1', 'A334'])

#Should I keep the protein regid?
#protein.flexibility['G649'].setLibraryRotamers(osprey.WILD_TYPE, 'TYR', 'ALA', 'VAL', 'ILE', 'LEU').addWildTypeRotamers().setContinuous()
#protein.flexibility['G650'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
#protein.flexibility['G651'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
#protein.flexibility['G654'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

## MFS
# * No mutable residues on the ligand
# * We allow translation and rotation of ligand
# * No translation rotation on protein
# * We set residues on the protein flexible

# 1F - 211Glu, 212Ser, 213Pro, 214Pro
# 2P - 205Tyr, 210Gly,
# 3L - 106Ile, 142Leu, 210Gly, 
# 4K - 179Glu
# 5R - 143Asp
# 6H - 177Thr, 178Pro, 179Glu
# 7D - 177Thr
# 8K - 139Asp, 160Leu
# 9V - 175Cys, 176Gly, 178Pro
# 10D - 174Phe, 175Cys
# 11D -  174Phe
# 12L - 174Phe, 178Pro, 181Leu, 186Leu, 223Phe
# 13S - 174Phe
#14K - 174Phe

#flex_residues = ['A211']

#for res in flex_residues:
    #protein.flexibility[res].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# Define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['B1', 'B14'])


#These residues, if changed, result in decreased binding capacity of MKI
flex_residues = ['B3', 'B5', 'B9', 'B14']

for res in flex_residues:
    ligand.flexibility[res].setLibraryRotamers(osprey.WILD_TYPE, 'GLY', 'ASN', 'ALA').addWildTypeRotamers().setContinuous()

# Add flexibility and potential mutations to key ligand residues
#ligand.flexibility['B2'].setLibraryRotamers(osprey.WILD_TYPE, 'PHE', 'TYR').addWildTypeRotamers().setContinuous()
#ligand.flexibility['B5'].setLibraryRotamers(osprey.WILD_TYPE, 'LEU', 'ILE', 'VAL').addWildTypeRotamers().setContinuous()
#ligand.flexibility['B8'].setLibraryRotamers(osprey.WILD_TYPE, 'ASP', 'GLU').addWildTypeRotamers().setContinuous()
#ligand.flexibility['B10'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
#ligand.flexibility['B12'].setLibraryRotamers(osprey.WILD_TYPE, 'ASN', 'GLN').addWildTypeRotamers().setContinuous()
#ligand.flexibility['B14'].setLibraryRotamers(osprey.WILD_TYPE, 'SER', 'THR').addWildTypeRotamers().setContinuous()

# make the conf space for the protein
proteinConfSpace = osprey.ConfSpace(protein)

# make the conf space for the ligand
ligandConfSpace = osprey.ConfSpace(ligand)

# make the conf space for the protein+ligand complex
complexConfSpace = osprey.ConfSpace([protein, ligand])

# Set flexiblity for the protein and ligand strands
trflex = osprey.c.confspace.StrandFlex.TranslateRotate(10, 2.0)

# Configure energy calculation
parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism)

# configure K*
kstar = osprey.KStar(
	proteinConfSpace,
	ligandConfSpace,
	complexConfSpace,
	epsilon=0.683,
	writeSequencesToFile=os.path.join(output_dir, 'mfs.kstar.results.tsv')
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

# Write K* scores to file
with open(os.path.join(output_dir, 'mfs.kstar.scores.tsv'), 'w') as f:
    f.write("sequence\tkstar_score\n")
    for scoredSequence in scoredSequences:
        f.write(f"{scoredSequence.sequence}\t{scoredSequence.score}\n")

# Save conformations for each scored sequence
for scoredSequence in scoredSequences:
    print("result:")
    print(f"\tsequence: {scoredSequence.sequence}")
    print(f"\tK* score: {scoredSequence.score}")

    numConfs = 10
    analysis = analyzer.analyze(scoredSequence.sequence, numConfs)
    print(analysis)

    pdb_filename = os.path.join(output_dir, f'seq.{scoredSequence.sequence}.pdb')
    analysis.writePdb(
        pdb_filename,
        f'Top {numConfs} conformations for sequence {scoredSequence.sequence}'
    )
    print(f"Saved top conformations to {pdb_filename}")
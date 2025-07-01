import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# read a PDB file for molecular info
mol = osprey.readPdb('/hpc/home/etm33/kinase_inhibitor_design/structures/processed/af3complex.clean.pdb') #pkn2.cleaned.filtered.pdb')

# define the protein strand-- here we want the entire protein(?)
protein = osprey.Strand(mol, templateLib=templateLib, residues=['A1', 'A334'])

## MFS
# * No mutable residues on the ligand
# * We allow translation and rotation of ligand
# * No translation rotation on protein
# * We set residues on the protein flexible

#1F - NA
#2P – GLY210
#3L - NA
#4K - NA
#5R – ASP143
#6H – NA
#7D – HIS24
#8K – PHE25
#9V – NA
#10D – PHE25, CYS175
#11D – CYS175
#12L – PHE174
#13S – NA
#14K – NA


flex_residues = ['A210', 
				 'A143', 
				 'A24', 
    			 'A25',
                 'A175',
                 'A174']

for res in flex_residues:
    protein.flexibility[res].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# Define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['B1', 'B14'])

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
	writeSequencesToFile='mfs.kstar.results.tsv'
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

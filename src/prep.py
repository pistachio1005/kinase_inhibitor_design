## prep.py, received from Henry Childs Jun 2025

import osprey
osprey.start()

# import the prep module after starting Osprey
import osprey.prep


# Before we get started:
# Preparing a molecule for Osprey can be a very error-prone process.
# Automated scripts, by their very nature, don't provide any feedback
# to the user between processing steps, so these errors are very likely
# to go unnoticed. For that reason, we strongly recommend that molecule
# preparations be done interactively using the GUI rather than an automated
# script. Nevertheless, if you already know where the errors in the prep
# process will happen, you can automate the prep using a script like this
# one, after customizing it to your specific situation.


# let's load a PDB file
# this one is a crystal structure of the Bovine PTPase complexed with HEPES
pdb_path = '8gaj.pdb'
pdb = osprey.prep.loadPDB(open(pdb_path, 'r').read())

# what molecules did we get?
print('Loaded %d molecules:' % len(pdb))
for mol in pdb:
    print('\t%s: %s' % (mol, osprey.prep.molTypes(mol)))

# looks like the PDB file has a protein chain and a small molecule
# the protein chain must be PTPase, and the small molecule must be HEPES
# (ignore the solvent molecules, if any)
target = pdb[0]
ligand = pdb[1]
mols = [target, ligand]

# start the local service that calls AmberTools for us
# NOTE: this will only work on Linux machines
with osprey.prep.LocalService():

    # Molecule Preparation Step 1: remove duplicate atoms
    # Duplicate atoms don't usually appear in files from the PDB,
    # but these errors can happen sometimes in modified PDB files.
    for mol in mols:
        # remove all but the first duplicated atom from each group
        for group in osprey.prep.duplicateAtoms(mol):
            for atomi in range(1, len(group.getAtoms())):
                group.remove(atomi)
                print('removed duplicate atom %s' % group)

    # Molecule Preparation Step 2: add missing heavy atoms
    # Somtimes atom positions are not well-resolved in the electron density,
    # or protein chain end-caps are missing.
    # But we still want to include these atoms in the molecular models to
    # be able to infer bonds correctly in the next step.
    for mol in mols:
        for missing_atom in osprey.prep.inferMissingAtoms(mol):
            missing_atom.add()
            print('added missing atom: %s' % missing_atom)

    # Molecule Preparation Step 3: add bonds
    # PDB files contain no explicit information about bonds, so we have to
    # infer where they might be based on the atoms we can see.
    for mol in mols:
        bonds = osprey.prep.inferBonds(mol)
        for bond in bonds:
            mol.getBonds().add(bond)
        print('added %d bonds to %s' % (len(bonds), mol))

    # Molecule Preparation Step 4: add hydrogens
    # Most PDB files based on X-ray crystallography don't describe Hydrogen atoms,
    # usually because light atoms aren't well-resolved in the electron density.
    # Generally, automatic protonation inference works well on proteins,
    # so let's just protonate PTPase for now, we'll protonate HEPES later
    for mol in mols:
        protonated_atoms = osprey.prep.inferProtonation(mol)
        for protonated_atom in protonated_atoms:
            protonated_atom.add()
        print('added %d hydrogens to %s' % (len(protonated_atoms), mol))

    # Moleclue Preparation Step 7: save the results
    pdb_path = '8gaj-ligand.pdb'
    open(pdb_path, 'w').write(osprey.prep.savePDB(ligand))
    print('saved prepared PDB to %s' % pdb_path)


print('Moleclue preparation complete!')

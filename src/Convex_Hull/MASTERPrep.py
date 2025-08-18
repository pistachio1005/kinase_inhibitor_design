from Bio.PDB import PDBParser, PDBIO
import sys
import warnings
import os
import fileinput

warnings.simplefilter('ignore')

def complex_renamer(complex_pdb: str, desired_target_id: str, desired_design_id: str):

    print("\nRenaming %s IDs. Target ID: %s, Design ID: %s.\n" % (complex_pdb, desired_target_id, desired_design_id))

    # get the old IDs
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", complex_pdb)
    model = structure[0]

    chain_ids = []
    for chain in model:
        chain_ids.append(chain.id)
    old_target_chain = chain_ids[0]
    old_peptide_chain = chain_ids[1]

    # change both IDs
    for chain in model:
        if chain.id == old_target_chain:
            chain.id = desired_target_id
        elif chain.id == old_peptide_chain:
            chain.id = desired_design_id

    # write out the new PDB with suffix -renamed.pdb
    io = PDBIO()
    io.set_structure(structure)
    filename = complex_pdb[:-4] + "-renamed.pdb"
    io.save(filename)
    print("Saved renamed PDB to %s" % filename)

    return filename

def reflect_design(complex_pdb: str, targetID: str, designID: str):
    print("Reflecting chain %s" % designID)

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", complex_pdb)
    model = structure[0]

    for chain in model:
        if chain.id == designID:
            counter = 1
            for res in chain:
                for a in res:
                    a.coord[2] = a.coord[2] * -1
                res.id = (' ', counter, ' ')
                counter += 1

    model.detach_child(targetID)
    io = PDBIO()
    io.set_structure(structure)
    filename = complex_pdb[:-4] + "-inverted-design.pdb"
    io.save(filename)
    print("Saved renamed PDB to %s" % filename)


def chain_relabel(input_directory: str, chainID: str):

    for pdb_path in os.listdir(input_directory):
        parser = PDBParser(PERMISSIVE=1)
        full_pdb_path = os.path.join(input_directory, pdb_path)
        M_pep_structure = parser.get_structure('match', full_pdb_path)
        M_peptide_model = M_pep_structure[0]

        for chain in M_peptide_model:
            chain.id = chainID

        io = PDBIO()
        io.set_structure(M_peptide_model)
        io.save(full_pdb_path)

        print("Renamed %s chain to %s" % (full_pdb_path, chainID))

# note: atom labels have not yet been corrected, so some side chain H are not correct
# we only need bb H atoms, so this is ok
def protonate_chains(matches_directory: str, target_pdb: str):

    import osprey
    osprey.start()
    import osprey.prep

    # protonate the target
    print("\n\n--Protonating L-target--\n\n")
    loaded_target_pdb = osprey.prep.loadPDB(open(target_pdb, 'r').read())
    targetChain = loaded_target_pdb[0]
    D_DesignChain = loaded_target_pdb[1]

    with osprey.prep.LocalService():
        osprey.prep.deprotonate(targetChain)
        protonated_atoms = osprey.prep.inferProtonation(targetChain)
        for protonated_atom in protonated_atoms:
            protonated_atom.add()
        print('Added %d hydrogens to %s' % (len(protonated_atoms), targetChain))

    # writing PDB via OSPREY is a pain, so we'll do some messy manual work to add chain sequentially
    open(target_pdb, 'w').write(osprey.prep.savePDB(targetChain))
    open(target_pdb, 'a').write(osprey.prep.savePDB(D_DesignChain))
    for line in fileinput.input(target_pdb, inplace=True):
        if "END" in line:
            print("TER")
        elif "REMARK" in line:
            continue
        else:
            print(line, end='')
    print("saved protonated PDB to %s" % target_pdb)

    # protonate each match (in-place)
    print("\n\n--Protonating L-matches--\n\n")
    for pdb_path in os.listdir(matches_directory):
        full_pdb_path = os.path.join(matches_directory, pdb_path)
        prep_pdb = osprey.prep.loadPDB(open(full_pdb_path, 'r').read())

        designChain = prep_pdb[0]

        with osprey.prep.LocalService():
            osprey.prep.deprotonate(designChain)
            protonated_atoms = osprey.prep.inferProtonation(designChain)
            for protonated_atom in protonated_atoms:
                protonated_atom.add()
            print('Added %d hydrogens to %s' % (len(protonated_atoms), designChain))

        open(full_pdb_path, 'w').write(osprey.prep.savePDB(designChain))
        print("saved protonated PDB to %s" % full_pdb_path)


def scaffold_generator(complex_pdb: str, designID: str, input_directory: str, output_directory: str):

    print("\n\n------ running scaffold generator ------\n\n")

    print("Storing outputs in %s" % output_directory)
    os.mkdir(output_directory)

    total_disjoint = 0

    for pdb_path in os.listdir(input_directory):

        # get the L-L complex
        parser = PDBParser(PERMISSIVE=1)
        complex_structure = parser.get_structure('complex', complex_pdb)
        complex_model = complex_structure[0]

        # get the MASTER L-peptide
        full_pdb_path = os.path.join(input_directory, pdb_path)
        M_pep_structure = parser.get_structure('match', full_pdb_path)
        M_peptide_model = M_pep_structure[0]

        # skip the match if it has a disjoint segment
        has_disjoint = False
        last_res = 0
        for chain in M_peptide_model:
            for res in chain:
                if last_res == 0:
                    last_res = res.id[1]
                elif (res.id[1] - last_res) != 1:
                    has_disjoint = True
                    print("SKIP: %s has a disjoint segment and will not be prepared" % pdb_path)
                    total_disjoint += 1
                    break
                else:
                    last_res = res.id[1]
        if has_disjoint:
            continue

        # flip each MASTER L-match to D-space, and change numbering and ID (if needed)
        for chain in M_peptide_model:
            counter = 1
            for res in chain:
                for a in res:
                    a.coord[2] = a.coord[2] * -1
                res.id = (' ', counter, ' ')
                counter += 1

        # MASTER already aligns, so just add the D-pep
        for chain in M_peptide_model:
            complex_model.add(chain)

        # remove the old L-peptide from the complex
        complex_model.detach_child(designID)

        # print new PDB to output_directory
        io = PDBIO()
        io.set_structure(complex_model)
        pdb_name = pdb_path[:-4] + "-complex.pdb"
        savename = os.path.join(output_directory, pdb_name)
        io.save(savename)

        # update terminal
        print("------ completed %s ------" % pdb_path)

    print("\n\nscaffold generation completed")
    print("Deleted %s matches due to disjoint segments\n\n" % total_disjoint)

## KStarPrep.py, received from Henry Childs Jun 2025

import math
import time

from Bio.PDB import PDBParser, PDBIO
from scipy.spatial import ConvexHull as ScipyConvexHull
from scipy.spatial import Delaunay

import os
import shutil

# BioPython throws errors that aren't really errors, so we'll ignore
import warnings

warnings.simplefilter('ignore')

from Find_Doublets import singlechain_design_info
from Confspace_Combiner import combine_Confspaces

# startup OSPREY
import osprey

osprey.start()
import osprey.prep


# class that holds confspace info specs
class ConfSpaceSpecs:
    def __init__(self, doublet, flexset, mutations, graph_path):
        self.doublet = doublet
        self.flexset = flexset
        self.mutations = mutations
        self.graph_path = graph_path


def doublet_confspace_info(pdb_filename, designID, design_muts: list, design_chirality, outHulls, FixedResidues: list):
    print("--Obtaining confspace info using hulls--")

    doublet, flex = singlechain_design_info(pdb_filename, outHulls, designID, design_muts, True, design_chirality, FixedResidues)

    all_pairs = []
    all_confspace = []

    for d in doublet:
        all_pairs.append(d)

    for p in all_pairs:
        if len(p) == 2:
            doub = list(p)
            res1_flex = flex[doub[0] - 1]
            res2_flex = flex[doub[1] - 1]
            flex_set = list(set(res1_flex).union(res2_flex))
            newspace = ConfSpaceSpecs(doub, flex_set, design_muts, doublet)
            all_confspace.append(newspace)
        elif len(p) == 1:
            # todo make sure this works
            doub = list(p)
            res_flex = flex[doub[0] - 1]
            newspace = ConfSpaceSpecs(doub, res_flex, design_muts, doublet)
            print("CHECKME! Found an island")
            exit()
        else:
            print("ERROR! 0 or > 2 residues in doublet")

    return all_confspace


def array_hull_coords(pdb_name: str):

    xyz_coords = []

    f = open(pdb_name, "r")
    for line in f:
        if "TER" in line:
            continue
        coords = line.split('     ')[2]
        coords_format = ""
        counter = -1
        for i in coords:
            if i == ' ' and counter == 0:
                coords_format += ","
                counter += 1
            if i in ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ".", "-"]:
                counter = 0
                coords_format += i

        all_coords = coords_format.split(',')
        atom_coords = []
        for i in all_coords:
            atom_coords.append(float(i))

        xyz_coords.append(atom_coords)

    return xyz_coords


def reduce_doublets(pdb_name: str, designID: str, targetID: str, specs: list, max_flex: int, method: str, hull_folder: str):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('complex', pdb_name)
    model = structure[0]

    trimmed_specs = []

    for doub in specs:
        if len(doub.flexset) > max_flex:
            print("Doublet %s has %s flexible residues which is > max flex of %s residues." % (doub.doublet, len(doub.flexset), max_flex))

            if method == "backbone_distance":
                print("Reducing using distance of CA in backbone")

                res1_coords = get_CA_coords(designID, doub.doublet[0], pdb_name)
                res2_coords = get_CA_coords(designID, doub.doublet[1], pdb_name)

                targetID = [c.id for c in model.get_chains() if c.id != designID]
                target_distances = dict()

                for t in doub.flexset:
                    target_coords = get_CA_coords(targetID[0], t, pdb_name)
                    res1_diff = euclidean_dist(res1_coords, target_coords)
                    res2_diff = euclidean_dist(res2_coords, target_coords)
                    if res1_diff < res2_diff:
                        target_distances[t] = res1_diff
                    else:
                        target_distances[t] = res2_diff

                ordered_targets = {k: v for k, v in sorted(target_distances.items(), key=lambda item: item[1])}
                print("CA distances for doublet %s:" % doub.doublet)
                print(ordered_targets)

                new_flex_set = []
                counter = 1
                for key in ordered_targets:
                    if counter > max_flex:
                        break
                    else:
                        new_flex_set.append(key)
                        counter += 1

                all_mutations = doub.mutations
                doub_path = doub.graph_path
                new_doublet = ConfSpaceSpecs(doub.doublet, new_flex_set, all_mutations, doub_path)
                trimmed_specs.append(new_doublet)
                print("Added reduced doublet %s with flex set %s" % (new_doublet.doublet, new_doublet.flexset))

            elif method == "volume_difference":
                print("Reducing using volume increase of composite hull")

                res1id = doub.doublet[0]
                res2id = doub.doublet[1]

                res1_file = ("%s/Chain%sRes%s.pdb" % (hull_folder, designID, res1id))
                res2_file = ("%s/Chain%sRes%s.pdb" % (hull_folder, designID, res2id))

                res1_coords = array_hull_coords(res1_file)
                res2_coords = array_hull_coords(res2_file)

                target_volumes = {}

                for target_resid in doub.flexset:
                    largest_volume = ["", float('inf')]
                    tar_filename = ("%s/Chain%sRes%s.pdb" % (hull_folder, targetID, target_resid))
                    target_coords = array_hull_coords(tar_filename)
                    hull1_mesh = Delaunay(res1_coords)
                    hull2_mesh = Delaunay(res2_coords)
                    for point in target_coords:
                        if hull1_mesh.find_simplex(point) >= 0:
                            res1_hull = ScipyConvexHull(res1_coords)
                            combined_hull = ScipyConvexHull(res1_coords + target_coords)
                            volume_diff = combined_hull.volume - res1_hull.volume
                            if volume_diff < largest_volume[1]:
                                largest_volume[0] = target_resid
                                largest_volume[1] = volume_diff

                        if hull2_mesh.find_simplex(point) >= 0:
                            res2_hull = ScipyConvexHull(res2_coords)
                            combined_hull = ScipyConvexHull(res2_coords + target_coords)
                            volume_diff = combined_hull.volume - res2_hull.volume
                            if volume_diff < largest_volume[1]:
                                largest_volume[0] = target_resid
                                largest_volume[1] = volume_diff

                    if largest_volume[0] == "":
                        print("ERROR! Something went wrong with flexible set reduction")
                        exit()

                    target_volumes[largest_volume[0]] = largest_volume[1]

                ordered_targets = {k: v for k, v in sorted(target_volumes.items(), key=lambda item: item[1])}
                print("Volume increases are: %s" % ordered_targets)

                new_flex_set = []
                counter = 1
                for key in ordered_targets:
                    if counter > max_flex:
                        break
                    else:
                        new_flex_set.append(key)
                        counter += 1

                all_mutations = doub.mutations
                doub_path = doub.graph_path
                new_doublet = ConfSpaceSpecs(doub.doublet, new_flex_set, all_mutations, doub_path)
                trimmed_specs.append(new_doublet)
                print("Added reduced doublet %s with flex set %s" % (new_doublet.doublet, new_doublet.flexset))

            elif method == "area_difference":
                print("Reducing using surface area change")

                res1id = doub.doublet[0]
                res2id = doub.doublet[1]

                res1_file = ("%s/Chain%sRes%s.pdb" % (hull_folder, designID, res1id))
                res2_file = ("%s/Chain%sRes%s.pdb" % (hull_folder, designID, res2id))

                res1_coords = array_hull_coords(res1_file)
                res2_coords = array_hull_coords(res2_file)

                res1_hull = ScipyConvexHull(res1_coords)
                res2_hull = ScipyConvexHull(res2_coords)

                res1_area = res1_hull.area
                res2_area = res2_hull.area

                res1_mesh = Delaunay(res1_coords)
                res2_mesh = Delaunay(res2_coords)

                target_volumes = {}

                for target_resid in doub.flexset:
                    largest_volume = ["", 0]
                    tar_filename = ("%s/Chain%sRes%s.pdb" % (hull_folder, targetID, target_resid))
                    target_coords = array_hull_coords(tar_filename)
                    target_hull = ScipyConvexHull(target_coords)
                    target_area = target_hull.area
                    for point in target_coords:
                        if res1_mesh.find_simplex(point) >= 0:
                            composite_hull = ScipyConvexHull(res1_coords + target_coords)
                            composite_area = composite_hull.area
                            area_diff = ((target_area + res1_area) - composite_area) / target_area
                            if area_diff > largest_volume[1]:
                                largest_volume[0] = target_resid
                                largest_volume[1] = area_diff
                        if res2_mesh.find_simplex(point) >= 0:
                            composite_hull = ScipyConvexHull(res2_coords + target_coords)
                            composite_area = composite_hull.area
                            area_diff = ((target_area + res2_area) - composite_area) / target_area
                            if area_diff > largest_volume[1]:
                                largest_volume[0] = target_resid
                                largest_volume[1] = area_diff

                    # small overlaps sometimes go negative, so just add as 0
                    if largest_volume[0] == "":
                        largest_volume[0] = target_resid
                        largest_volume[1] = 0

                    target_volumes[largest_volume[0]] = largest_volume[1]

                ordered_targets = {k: v for k, v in reversed(sorted(target_volumes.items(), key=lambda item: item[1]))}
                print("Area reductions are: %s" % ordered_targets)

                new_flex_set = []
                counter = 1
                for key in ordered_targets:
                    if counter > max_flex:
                        break
                    else:
                        new_flex_set.append(key)
                        counter += 1

                all_mutations = doub.mutations
                doub_path = doub.graph_path
                new_doublet = ConfSpaceSpecs(doub.doublet, new_flex_set, all_mutations, doub_path)
                trimmed_specs.append(new_doublet)
                print("Added reduced doublet %s with flex set %s" % (new_doublet.doublet, new_doublet.flexset))

        else:
            print("Doublet %s only has %s flexible residues, so will not be reduced" % (doub.doublet, len(doub.flexset)))
            trimmed_specs.append(doub)

    return trimmed_specs


def euclidean_dist(coord1, coord2):

    x_diff = (coord1[0] - coord2[0]) ** 2
    y_diff = (coord1[1] - coord2[1]) ** 2
    z_diff = (coord1[2] - coord2[2]) ** 2

    return math.sqrt((x_diff + y_diff + z_diff))


def get_CA_coords(chainID: str, resID: int, pdb_filename: str):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('complex', pdb_filename)
    model = structure[0]

    for chain in model:
        if chain.id == chainID:
            for residue in chain:
                if residue.id[1] == resID:
                    for atom in residue:
                        if atom.name == 'CA':
                            coords = atom.get_vector()
                            return list(coords)

    print("ERROR! Couldn't find CA for residue %s in chain %s" % (resID, chainID))
    exit()

def error_check_pdb(pdb_name: str):
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('complex', pdb_name)
    model = structure[0]

    # return error if N/C-term PRO present (OSPREY anchor chords can't handle these currently)
    for chain in model:
        chain_length = len(chain)
        counter = 0
        for res in chain:
            counter += 1
            if counter == 1 or counter == chain_length:
                atype = res.get_resname()
                if atype == "PRO":
                    print("ERROR! PDB %s contains an N or C-term proline, which can't be used in OSPREY" % pdb_name)
                    print("Halting process for this file")
                    return "ERROR_TERM_PRO"

    # change CD -> CD1 labelling, which is common mislabel in PDBs for these residues
    for chain in model:
        for res in chain:
            resname = res.get_resname()
            resnum = res.id[1]
            if resname in ("ILE", "LEU", "PHE", "TRP", "TYR"):
                for a in res:
                    if a.fullname == ' CD ':
                        a.fullname = ' CD1'
                        print("Mislabel found with carbon label. Changing residue %s%s to CD1." % (resname, resnum))

    new_pdb_name = pdb_name.split(".")[0] + "-corrected.pdb"

    io = PDBIO()
    io.set_structure(structure)
    io.save(new_pdb_name)
    print("Now using corrected PDB file %s for future fileprep" % new_pdb_name)

    return new_pdb_name


def combine_PDBs(pdb_list, outfile_name, delete_old: bool, designID: str, flip_design: bool):
    with open(outfile_name, 'w') as outfile:
        for fname in pdb_list:
            with open(fname) as infile:
                for line in infile:
                    if "REMARK" in line:
                        pass
                    elif "END\n" in line:
                        outfile.write("TER\n")
                    else:
                        outfile.write(line)

    if delete_old:
        for pdb_file in pdb_list:
            os.remove(pdb_file)

    if flip_design:
        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure('complex', outfile_name)
        model = structure[0]
        d_pep = model[designID]
        for r in d_pep:
            for a in r:
                a.coord[2] = a.coord[2] * -1

        io = PDBIO()
        io.set_structure(structure)
        io.save(outfile_name)
        print("Flipped chain back to D space and saved to %s" % outfile_name)


def osprey_fileprep_preprocess(pdb_name, designID: str, chirality: str, mindesign: bool, add_protons: bool):

    changed_pdb = False

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('complex', pdb_name)
    model = structure[0]

    chain_ids = []
    for chain in model:
        chain_ids.append(chain.id)

    if chirality == 'D':
        print("Have D-space chain, so created PDB in L-space for fileprep")
        d_pep = model[designID]
        for r in d_pep:
            for a in r:
                a.coord[2] = a.coord[2] * -1

        io = PDBIO()
        io.set_structure(structure)
        L_filename = pdb_name.split(".")[0] + "-flipped-design.pdb"
        io.save(L_filename)
        print("Saved flipped chain to %s" % L_filename)
        pdb_name = L_filename
        changed_pdb = True

    prep_pdb = osprey.prep.loadPDB(open(pdb_name, 'r').read())
    if designID == chain_ids[1]:
        target = prep_pdb[0]
        design = prep_pdb[1]
    elif designID == chain_ids[0]:
        target = prep_pdb[1]
        design = prep_pdb[0]
    else:
        print("ERROR! Was unable to find chain %s in PDB %s" % (designID, pdb_name))

    mols = [target, design]

    with osprey.prep.LocalService():
        for mol in mols:
            for group in osprey.prep.duplicateAtoms(mol):
                for atomi in range(1, len(group.getAtoms())):
                    group.remove(atomi)
                    changed_pdb = True
                    print('Removed duplicate atom %s' % group)

        for mol in mols:
            for missing_atom in osprey.prep.inferMissingAtoms(mol):
                missing_atom.add()
                changed_pdb = True
                print('Added missing atom: %s' % missing_atom)

        if add_protons:
            for mol in mols:
                osprey.prep.deprotonate(mol)
                protonated_atoms = osprey.prep.inferProtonation(mol)
                for protonated_atom in protonated_atoms:
                    protonated_atom.add()
                    changed_pdb = True
                print('Added %d hydrogens to %s' % (len(protonated_atoms), mol))

    if changed_pdb:
        new_pdb_design = pdb_name.split(".")[0] + "-design-temp.pdb"
        new_pdb_target = pdb_name.split(".")[0] + "-target-temp.pdb"
        new_pdb_complex = pdb_name.split("-")[0] + "-processed.pdb"
        open(new_pdb_design, 'w').write(osprey.prep.savePDB(design))
        open(new_pdb_target, 'w').write(osprey.prep.savePDB(target))
        if chirality == 'D':
            combine_PDBs([new_pdb_target, new_pdb_design], new_pdb_complex, True, designID, True)
        elif chirality != 'D':
            combine_PDBs([new_pdb_target, new_pdb_design], new_pdb_complex, True, designID, False)
        print("Note: corrected or added atoms in PDB %s. Created new PDB %s for future fileprep." % (pdb_name, new_pdb_complex))
        return new_pdb_complex

    # TODO minimize the target chain if requested

    return pdb_name


def make_omol_file(pdb_name, isDesign: bool):
    processed_pdb = osprey.prep.loadPDB(open(pdb_name, 'r').read())

    target = processed_pdb[0]
    design = processed_pdb[1]
    mols = [target, design]

    with osprey.prep.LocalService():
        for mol in mols:
            bonds = osprey.prep.inferBonds(mol)
            for bond in bonds:
                mol.getBonds().add(bond)
            print('added %d bonds to %s' % (len(bonds), mol))

    if isDesign:
        omol_name = pdb_name.split(".")[0] + "-design.omol"
        open(omol_name, 'w').write(osprey.prep.saveOMOL([design]))
    else:
        omol_name = pdb_name.split(".")[0] + "-target.omol"
        open(omol_name, 'w').write(osprey.prep.saveOMOL([target]))

    print("Saved to %s" % omol_name)
    return omol_name


def get_omol_res_chain(omol_name: str, resnum: str):
    resType = ""
    chainID = ""

    search_name = ("[molecule.%s.polymer]" % 0)

    f = open(omol_name, 'r')
    found_start = False
    found_chain = False
    for line in f:
        if found_start:
            chainID = line[1]
            found_chain = True
            found_start = False
        elif found_chain:
            id_region = line.split(",")[0]
            curr_resnum = int(id_region.split("\"")[1])
            if curr_resnum == resnum:
                type_region = line.split(",")[1]
                curr_type = type_region.split("\"")[1]
                resType = curr_type
                break
        elif search_name in line:
            found_start = True

    return resType, chainID


def make_target_confspace(confspec: ConfSpaceSpecs, omol_name: str):
    osprey_omol = osprey.prep.loadOMOL(open(omol_name, 'r').read())
    target_conf = osprey_omol[0]
    target_conf_space = osprey.prep.ConfSpace(osprey_omol)
    lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'lovell2000-osprey3').load()
    target_conf_space.getConflibs().add(lovell2000)

    flex_res = confspec.flexset
    for resnum in flex_res:
        resType, chainID = get_omol_res_chain(omol_name, resnum)
        new_flex = target_conf_space.addPosition(osprey.prep.ProteinDesignPosition(target_conf, chainID, str(resnum)))
        target_conf_space.addMutations(new_flex, resType)

    print('Target confspace for doublet %s:' % confspec.doublet)
    for pos in target_conf_space.positions():
        print('\t%6s flexing: %s' % (pos.getName(), target_conf_space.getMutations(pos)))

    for pos in target_conf_space.positions():
        for mutation in target_conf_space.getMutations(pos):
            target_conf_space.addConformationsFromLibraries(pos, mutation)
        # if pos.getType() in target_conf_space.getMutations(pos):
        #     target_conf_space.addWildTypeConformation(pos)

    dihedral_settings = osprey.prep.DihedralAngleSettings()
    for pos in target_conf_space.positions():
        for mutation in target_conf_space.getMutations(pos):
            for conf_info in target_conf_space.getConformations(pos, mutation):
                for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                    conf_info.getMotions().add(motion)

    print(
        'Target conformation space describes %s conformations for doublet %s' % (target_conf_space.countConformations(),
                                                                                 confspec.doublet))

    curr_doublet = confspec.doublet
    doublet_print = ("[%s_%s]" % (curr_doublet[0], curr_doublet[1]))
    target_conf_path = omol_name.split(".")[0] + "-" + doublet_print + ".confspace"
    open(target_conf_path, 'w').write(osprey.prep.saveConfSpace(target_conf_space))
    print("Saved target confspace to %s" % target_conf_path)

    return target_conf_path

def get_omol_length(omol_name: str):

    search_name = ("[molecule.%s.polymer]" % 0)

    f = open(omol_name, 'r')
    at_start = False
    num_res = 0
    for line in f:
        if at_start:
            if "{" in line:
                num_res += 1
        if search_name in line:
            at_start = True

    if num_res == 0:
        print("ERROR! Unable to find chain length from OMOL file")
        exit()

    return num_res


def make_design_confspace(confspec: ConfSpaceSpecs, omol_name: str, mutTypes: list, translate_rotate: bool,
                          chirality: str, otherResAla: bool, fixed_residues: list):
    osprey_omol = osprey.prep.loadOMOL(open(omol_name, 'r').read())
    design_conf = osprey_omol[0]
    design_conf_space = osprey.prep.ConfSpace(osprey_omol)

    lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'D-lovell2000-osprey3').load()
    design_conf_space.getConflibs().add(lovell2000)

    preHis_muts = len(mutTypes)
    try:
        mutTypes.remove('HIP')
        mutTypes.remove('HID')
        mutTypes.remove('HIE')
    except:
        pass
    postHis_muts = len(mutTypes)
    if preHis_muts != postHis_muts:
        mutTypes.append("HIS")

    for resnum in confspec.doublet:
        design_length = get_omol_length(omol_name)
        resType, chainID = get_omol_res_chain(omol_name, resnum)
        new_flex = design_conf_space.addPosition(osprey.prep.ProteinDesignPosition(design_conf, chainID, str(resnum)))
        if resnum in fixed_residues:
            print("Design chain %s%s has a fixed identity, so only getting WT flexibility" % (resType, resnum))
            design_conf_space.addMutations(new_flex, resType)
        elif resnum != 1 and resnum != design_length:
            design_conf_space.addMutations(new_flex, mutTypes)
        else:
            print("N or C terminus in doublet, so excluding PRO from mutations (if requested)")
            no_PRO = []
            for a in mutTypes:
                if a != 'PRO':
                    no_PRO.append(a)
            design_conf_space.addMutations(new_flex, no_PRO)

    if otherResAla:
        print("Now changing non-mutant design chain residues to ALA")
        design_length = get_omol_length(omol_name)
        non_muts = []
        for i in range(1, design_length+1):
            if i != confspec.doublet[0] and i != confspec.doublet[1]:
                non_muts.append(i)
        for n in non_muts:
            resType, chainID = get_omol_res_chain(omol_name, confspec.doublet[0])
            new_ala = design_conf_space.addPosition(osprey.prep.ProteinDesignPosition(design_conf, chainID, str(n)))
            design_conf_space.addMutations(new_ala, 'ALA')

    print('Design confspace for doublet %s:' % confspec.doublet)
    for pos in design_conf_space.positions():
        print('\t%6s mutating: %s' % (pos.getName(), design_conf_space.getMutations(pos)))

    for pos in design_conf_space.positions():
        for mutation in design_conf_space.getMutations(pos):
            design_conf_space.addConformationsFromLibraries(pos, mutation)
        # if pos.getType() in design_conf_space.getMutations(pos):
        #     design_conf_space.addWildTypeConformation(pos)

    dihedral_settings = osprey.prep.DihedralAngleSettings()
    for pos in design_conf_space.positions():
        for mutation in design_conf_space.getMutations(pos):
            for conf_info in design_conf_space.getConformations(pos, mutation):
                for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                    conf_info.getMotions().add(motion)

    print("WARNING: translation/rotation currently on works with K* searches using JDK19 (see nom branch). Python "
          "versions (from main) will not run.")
    if translate_rotate:
        design_conf_space.addMotion(osprey.prep.moleculeTranslationRotation(design_conf))

    print(
        'Design conformation space describes %s conformations for doublet %s' % (design_conf_space.countConformations(),
                                                                                 confspec.doublet))

    curr_doublet = confspec.doublet
    doublet_print = ("[%s_%s]" % (curr_doublet[0], curr_doublet[1]))
    design_conf_path = omol_name.split(".")[0] + "-" + doublet_print + ".confspace"
    open(design_conf_path, 'w').write(osprey.prep.saveConfSpace(design_conf_space))
    print("Saved design confspace to %s" % design_conf_path)

    return design_conf_path


def compile_confspaces(spaces: list):
    for s in spaces:
        confspace = osprey.prep.loadConfSpace(open(s, 'r').read())
        save_path = s.split(".")[0] + ".ccsx"

        compiler = osprey.prep.ConfSpaceCompiler(confspace)

        compiler.getForcefields().add(osprey.prep.Forcefield.Amber96)
        compiler.getForcefields().add(osprey.prep.Forcefield.EEF1)

        print('Compiling %s' % s)
        progress = compiler.compile()
        progress.printUntilFinish(10000)
        report = progress.getReport()

        if report.getError() is not None:
            raise Exception('Compilation failed', report.getError())

        open(save_path, 'wb').write(osprey.prep.saveCompiledConfSpace(report.getCompiled()))
        print('Saved compiled confspace to %s' % save_path)


def organize_kstar_files(out_directory: str, doublet: ConfSpaceSpecs, target_name: str, design_name: str, complex_name: str, KStarBash: str):
    doublets = doublet.doublet
    doublet_name = ("[%s_%s]" % (doublets[0], doublets[1]))
    new_directory = ("%s/kstar-%s" % (out_directory, doublet_name))
    try:
        os.mkdir(new_directory)
    except Exception as e:
        print("ERROR! Unable to make directory %s" % new_directory)
        print(e)
        exit()

    t_dest = ("%s/%s" % (new_directory, target_name))
    print("Moving %s to %s" % (target_name, t_dest))
    os.rename(target_name, t_dest)

    d_dest = ("%s/%s" % (new_directory, design_name))
    print("Moving %s to %s" % (design_name, d_dest))
    os.rename(design_name, d_dest)

    c_dest = ("%s/%s" % (new_directory, complex_name))
    print("Moving %s to %s" % (complex_name, c_dest))
    os.rename(complex_name, c_dest)

    target_ccsx = target_name.split(".")[0] + ".ccsx"
    target_ccsx_path = ("%s/%s" % (new_directory, "target.ccsx"))
    print("Moving %s to %s" % (target_ccsx, target_ccsx_path))
    os.rename(target_ccsx, target_ccsx_path)

    design_ccsx = design_name.split(".")[0] + ".ccsx"
    design_ccsx_path = ("%s/%s" % (new_directory, "design.ccsx"))
    print("Moving %s to %s" % (design_ccsx, design_ccsx_path))
    os.rename(design_ccsx, design_ccsx_path)

    complex_ccsx = complex_name.split(".")[0] + ".ccsx"
    complex_ccsx_path = ("%s/%s" % (new_directory, "complex.ccsx"))
    print("Moving %s to %s" % (complex_ccsx, complex_ccsx_path))
    os.rename(complex_ccsx, complex_ccsx_path)

    try:
        e_direc = ("%s/%s" % (new_directory, "ensembles"))
        os.mkdir(e_direc)
    except Exception as e:
        print("ERROR! Unable to make directory %s" % e_direc)
        print(e)
        exit()

    shutil.copy(KStarBash, new_directory)


# function for file prepping D:L designs
def osprey_fileprep_kstar(pdb_input, out_directory: str, design_muts: list, confspecs: list, designID: str,
                          design_chirality: str,
                          minimize_design: bool,
                          add_protons: bool, translation_rotation: bool, otherResAla: bool, KStarBash: str, FixedResidues: list):
    print("Creating OSPREY files for your design")

    print("--Checking input PDB file for any labelling or content issues--")
    correct_pdb = error_check_pdb(pdb_input)
    if correct_pdb == "ERROR_TERM_PRO":
        return

    print("--Preprocessing PDB file--")
    preprocess_pdb = osprey_fileprep_preprocess(correct_pdb, designID, design_chirality, minimize_design, add_protons)

    print("--Making OMOL files--")
    design_omol = make_omol_file(preprocess_pdb, True)
    target_omol = make_omol_file(preprocess_pdb, False)

    print("--Making OSPREY confspace files--")
    for spec in confspecs:

        print("Now making confspace files for doublet %s with flexset %s" % (spec.doublet, spec.flexset))
        if design_chirality == 'L':
            print("TODO chirality is L for design")
            exit()
        elif design_chirality == 'D':
            target_conf_name = make_target_confspace(spec, target_omol)
            design_conf_name = make_design_confspace(spec, design_omol, design_muts, translation_rotation, design_chirality,
                                                 otherResAla, FixedResidues)
            complex_conf_name = combine_Confspaces(target_conf_name, design_conf_name, True)

        print("--Compiling OSPREY confspace files for doublet %s--" % spec.doublet)
        with osprey.prep.LocalService():
            compile_confspaces([target_conf_name, design_conf_name, complex_conf_name])

        print("--Organizing files for %s into directory %s--" % (spec.doublet, out_directory))
        organize_kstar_files(out_directory, spec, target_conf_name, design_conf_name, complex_conf_name, KStarBash)

    design_omol_path = ("%s/%s" % (out_directory, design_omol))
    target_omol_path = ("%s/%s" % (out_directory, target_omol))
    os.rename(design_omol, design_omol_path)
    os.rename(target_omol, target_omol_path)


def save_doublets(specs: list, filename: str):

    all_doublets = ""

    for i in range(0, len(specs)):
        if i != (len(specs) - 1):
            all_doublets += str(specs[i].doublet)
            all_doublets += ','
        else:
            all_doublets += str(specs[i].doublet)

    with open(filename, "w") as f:
        f.write(all_doublets)
        f.write('\n')
        for s in specs:
            new = ""
            new += str(s.doublet)
            new += ' '
            new += str(s.flexset)
            f.write(new)
            f.write('\n')

    print("Saved doublet info to %s" % filename)

    f.close()


def make_target_ABAA(omol_name: str, flexibility: list):
    osprey_omol = osprey.prep.loadOMOL(open(omol_name, 'r').read())
    target_conf = osprey_omol[0]
    target_conf_space = osprey.prep.ConfSpace(osprey_omol)
    lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'lovell2000-osprey3').load()
    target_conf_space.getConflibs().add(lovell2000)

    for resnum in flexibility:
        resType, chainID = get_omol_res_chain(omol_name, resnum)
        new_flex = target_conf_space.addPosition(osprey.prep.ProteinDesignPosition(target_conf, chainID, str(resnum)))
        target_conf_space.addMutations(new_flex, resType)

    print('Target confspace for ABAA:')
    for pos in target_conf_space.positions():
        print('\t%6s flexing: %s' % (pos.getName(), target_conf_space.getMutations(pos)))

    for pos in target_conf_space.positions():
        for mutation in target_conf_space.getMutations(pos):
            target_conf_space.addConformationsFromLibraries(pos, mutation)

    dihedral_settings = osprey.prep.DihedralAngleSettings()
    for pos in target_conf_space.positions():
        for mutation in target_conf_space.getMutations(pos):
            for conf_info in target_conf_space.getConformations(pos, mutation):
                for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                    conf_info.getMotions().add(motion)

    print(
        'Target conformation space describes %s conformations' % (target_conf_space.countConformations()))

    target_conf_path = omol_name.split(".")[0] + "-[ABAA].confspace"
    open(target_conf_path, 'w').write(osprey.prep.saveConfSpace(target_conf_space))
    print("Saved target confspace to %s" % target_conf_path)

    return target_conf_path


def target_flexibility_ABAA(pdb_name: str, designID: str, design_chirality: str):
    print("Finding flexibility for ABAA:")
    all_specs = singlechain_design_info(pdb_name, "", designID, ['VAL'], False, design_chirality, [])

    target_flex = []

    for flexlist in all_specs[1]:
        for res in flexlist:
            if res not in target_flex:
                target_flex.append(res)

    print("Flexible target residues are: %s" % target_flex)

    return target_flex


def find_gly_pro_index(pdb_name: str, chain_id: str):

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure('complex', pdb_name)
    model = structure[0]

    gly_pro_index = []

    for chain in model:
        if chain.id == chain_id:
            for residue in chain:
                resname = residue.resname
                resid = residue.id[1]
                if resname == 'GLY':
                    gly_pro_index.append(resid)
                elif resname == 'PRO':
                    gly_pro_index.append(resid)

    print("Not mutating residues %s to ALA because GLY or PRO" % gly_pro_index)

    return gly_pro_index


def make_design_ABAA(omol_name: str, designID: str, pdb_name: str):
    osprey_omol = osprey.prep.loadOMOL(open(omol_name, 'r').read())
    design_conf = osprey_omol[0]
    design_conf_space = osprey.prep.ConfSpace(osprey_omol)

    lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'D-lovell2000-osprey3').load()
    design_conf_space.getConflibs().add(lovell2000)

    design_length = get_omol_length(omol_name)
    gly_pro_residues = find_gly_pro_index(pdb_name, designID)
    for i in range(1, design_length+1):
        if i not in gly_pro_residues:
            new_ala = design_conf_space.addPosition(osprey.prep.ProteinDesignPosition(design_conf, designID, str(i)))
            design_conf_space.addMutations(new_ala, 'ALA')

    print('Design confspace for design:')
    for pos in design_conf_space.positions():
        print('\t%6s mutating: %s' % (pos.getName(), design_conf_space.getMutations(pos)))

    for pos in design_conf_space.positions():
        for mutation in design_conf_space.getMutations(pos):
            design_conf_space.addConformationsFromLibraries(pos, mutation)

    dihedral_settings = osprey.prep.DihedralAngleSettings()
    for pos in design_conf_space.positions():
        for mutation in design_conf_space.getMutations(pos):
            for conf_info in design_conf_space.getConformations(pos, mutation):
                for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                    conf_info.getMotions().add(motion)

    print("WARNING: translation/rotation currently on works with K* searches using JDK19 (see nom branch). Python "
          "versions (from main) will not run.")
    design_conf_space.addMotion(osprey.prep.moleculeTranslationRotation(design_conf))

    print(
        'Design conformation space describes %s conformations' % (design_conf_space.countConformations()))

    design_conf_path = omol_name.split(".")[0] + "-[ABAA].confspace"
    open(design_conf_path, 'w').write(osprey.prep.saveConfSpace(design_conf_space))
    print("Saved design confspace to %s" % design_conf_path)

    return design_conf_path


def osprey_fileprep_ABAA(pdb_input: str, out_directory: str, designID: str, design_chirality: str, minimize_design: bool, add_protons: bool, KStarBash: str):

    print("Creating OSPREY files for ABAA")

    print("--Checking input PDB file for any labelling or content issues--")
    correct_pdb = error_check_pdb(pdb_input)
    if correct_pdb == "ERROR_TERM_PRO":
        shutil.rmtree(out_directory)
        return

    print("--Preprocessing PDB file--")
    preprocess_pdb = osprey_fileprep_preprocess(correct_pdb, designID, design_chirality, minimize_design, add_protons)

    print("--Making OMOL files--")
    design_omol = make_omol_file(preprocess_pdb, True)
    target_omol = make_omol_file(preprocess_pdb, False)

    print("--Determining flexible set--")
    target_flex = target_flexibility_ABAA(pdb_input, designID, design_chirality)

    print("--Making confspace files--")
    target_conf_name = make_target_ABAA(target_omol, target_flex)
    design_conf_name = make_design_ABAA(design_omol, designID, preprocess_pdb)
    complex_conf_name = combine_Confspaces(target_conf_name, design_conf_name, True)

    print("--Compiling OSPREY confspace files")
    with osprey.prep.LocalService():
        compile_confspaces([target_conf_name, design_conf_name, complex_conf_name])

    print("--Organizing files into %s--" % out_directory)
    bbonly_name = ("%s-bbonly.pdb" % (pdb_input.split(".")[0]))
    os.rename(bbonly_name, ("%s/%s" % (out_directory, bbonly_name)))


    os.rename(correct_pdb, ("%s/%s" % (out_directory, correct_pdb)))
    os.rename(preprocess_pdb, ("%s/%s" % (out_directory, preprocess_pdb)))
    flipped_name = ("%s-flipped-design.pdb" % (correct_pdb.split(".")[0]))

    try:
        os.rename(flipped_name, ("%s/%s" % (out_directory, flipped_name)))
    except Exception:
        pass

    os.rename(design_omol, ("%s/%s" % (out_directory, design_omol)))
    os.rename(target_omol, ("%s/%s" % (out_directory, target_omol)))

    os.rename(target_conf_name, ("%s/%s" % (out_directory, target_conf_name)))
    os.rename(design_conf_name, ("%s/%s" % (out_directory, design_conf_name)))
    os.rename(complex_conf_name, ("%s/%s" % (out_directory, complex_conf_name)))

    os.rename(target_conf_name.split('.')[0] + ".ccsx", ("%s/%s" % (out_directory, "target.ccsx")))
    os.rename(design_conf_name.split('.')[0] + ".ccsx", ("%s/%s" % (out_directory, "design.ccsx")))
    os.rename(complex_conf_name.split('.')[0] + ".ccsx", ("%s/%s" % (out_directory, "complex.ccsx")))

    os.mkdir(("%s/%s" % (out_directory, "ensembles")))
    shutil.copy(KStarBash, out_directory)

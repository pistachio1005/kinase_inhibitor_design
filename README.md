# Kinase Inhibitor Design

## Current Workflow

### Molecule Preparation

This project requires models in `.pdb` format, and any modifications to the models should be done _before_ using the scripts in this repo.

* `sbatch prep.sh [path_to_pdb]`

This script will identify missing atoms/residues, add them, and then hydrogenate the protein (Assuming pH=7.0).

Note: If phosphorylated amino acids in the protein structure, `phosphorylated_to_canonical.py` will convert those residues to their canonical form.

Note that we have the follwoing directories for storing structures:

```
kinase_inhibitor_design
  src/
     /structures/[normal structures go here]
                 processed/[processed structures go here]
```

### Minimal Flexible Set ("Fitting" the Ligand to Protein)

* `sbatch MFS_SCOPE.sh [path_to_pdb]`

The current MFS framework sets flexibility based on `SCOPE` using statistics to pair down the number of residues we would set to flexible.

The current method:

* Sets **ligand** residues to flexible if:
    * They had WT rotamer orientations that could clash
    * The ligand residues clash with an “above average” number of residues on the protein
        * We take the average number of clashes with protein residues for each ligand residue to compute this

* Set **protein** residues to flexible if:
    * They have possible clashes with 2 or more residues on the ligand

### Positive Design

* `sbatch MFS_Design.sh [path_to_pdb]`

This is a work in progress.

## Installing OPSREY on Ubuntu (local machine) with IntelliJ

Using the terminal, in the `~/IdeaProjects` directory:

```
git clone https://github.com/donaldlab/OSPREY3.git
```

In idea, click the settings (gear) button, and then the "Project Structure" button to set the JDK to: 17 Oracle OpenJDK 17.0.12. Be sure to click "apply"

Note: If you look at "Build", OSPREY3 should begin building under these conditions

Once building is finished, click the Gradle Icon (on the right side of the IntelliJ IDE) and you will see the directories:

```
>Tasks
>Dependencies
>buildSrc
```

We run:

```
Tasks > build > pythonWheel
```

Result should be something along the lines of:

```
BUILD SUCCESSFUL in 2m 44s
```

## A Note on `Convex_Hull`

The `Convex_Hull` folder points to a (likely deprecated) version of the following github: [Convex_Hull](https://github.com/henry-childs/Convex_Hull).

In this project, we utilize `Convex_Hull`'s "Side Chain Orientation and Position Evaluation" (SCOPE) tool, written by [Henry Childs](https://github.com/henry-childs).

For this to be possible, we place the `Convex_Hull` directory into our `kinase_inhibitor_design/src` folder, which provides us with convenient access to `Convex_Hull`'s modules.

## Bonus: Conda Setup

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```

**Environment Creation/Activation**

```
conda create env -n Osprey3
conda activate Opsrey3
```

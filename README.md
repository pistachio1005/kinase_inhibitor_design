# kinase_inhibitor_design

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

Result should be:

```
BUILD SUCCESSFUL in 2m 44s
```

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
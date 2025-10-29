# ===============================================================
# Bond distances
# ===============================================================
# Finds bond distance between all atoms in a .xyz file
# and creates a dictionary with information of all atoms in 
# the molecule
# The dictionary is exported as a .json file
# Multiple .xyz files can be ran at the same time
#
# File should be ran from the same folder as the .xyz files
# New files will be created in the same folder as well, so if
# that is not desired, change the code accordingly 
# ===============================================================

import re
import json

files = ["ethanol", "sif4-pbe-7-6-opt"] # without .xyz
molecules = {}


def format_molecule_to_dict(file):
    """
    Input: name of .xyz file
    Output: nested dictionary of all atoms in the file
    Example output:
    {C1: {"element" : "C", "x" : 0.01247000, "y" = 0.02254000, "z" : 1.08262000},
     C2: {"element" : "C", "x" : 0.02024000, "y" = 0.02254000, "z" : 1.08262000}}
    """
    # List of lists where each list is on line from the .xyz file
    atoms_raw = []
    # List of dicts with formatted information
    atoms_collected = []
    # The dict that collects all labeled atoms
    mol = {}
    
    # Reads last geometry of xyz file and ads to atoms_raw
    with open(f"{file}.xyz", "r") as f:
        geomfile = f.readlines()
        n_atoms = int(geomfile[0])
        for i in range(1, n_atoms + 1):
            line = geomfile[-i].strip()
            atoms_raw.append(line)

    # Reformats lines in atoms_raw to a dictionary
    for atom in atoms_raw:
        mol_data = re.search(r"([A-Z][a-z]?)\s+(-?\d+\.?(?:\d+)?(?:e-?\d+)?)\s+(-?\d+\.?(?:\d+)?(?:e-?\d+)?)\s+(-?\d+\.?(?:\d+)?(?:e-?\d+)?)", atom)

        x, y, z = float(mol_data.group(2)), float(mol_data.group(3)), float(mol_data.group(4))
        element = mol_data.group(1)
    
        atoms_collected.append({"element": element, "x": x, "y": y, "z": z})

    # Reveres list as to get atomes ordered from top to match xyz file
    atoms_collected.reverse()

    # Counts number of each element, labels them and combines all atoms for the final dict
    counts = {}
    for atom in atoms_collected:
        element = atom["element"]
        counts[element] = counts.get(element, 0) + 1
        atom_name = f"{element}{counts[element]}"
        mol[atom_name] = atom
    return mol


def find_distance(atom1, atom2):
    """
    Input: two atoms
    Output: bond distance
    """
    x1 = atom1["x"]
    y1 = atom1["y"]
    z1 = atom1["z"]
    x2 = atom2["x"]
    y2 = atom2["y"]
    z2 = atom2["z"]
    return ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**0.5

def write_distance_file(d_file, mol_dict):
    """
    Input: name of file where distances will be written (without file extension)
    Output: none
    
    Writes a text file listing all bond distances
    of a given molecule dictionary
    New file will be called <name of .xyz file>_distances.txt
    """
    atoms = mol_dict
    atom_names = list(atoms.keys())

    
    with open(f"{d_file}_distances.txt", "w") as f:
        f.write(f"Distances for {d_file} (Same unit as .xyz file):\n\n")

        for i in range(len(atom_names)):
            for j in range(i + 1, len(atom_names)):
                atom1, atom2 = atom_names[i], atom_names[j]
                d = find_distance(atoms[atom1], atoms[atom2])
                f.write(f"{atom1} - {atom2}: {d}\n")


# Loops over lists of files and creates one distance filer per .xyz file
for file in files:
    molecules[file] = format_molecule_to_dict(file)
    write_distance_file(file, molecules[file])

# Creates json file containing dictionary with molecule information
with open("molecule_data.json", "w") as f:
    f.write(json.dumps(molecules, indent=4))

print("Done :)")
import math
import csv
from typing import List, Tuple, Dict
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
import os


def get_file_path(file_description):
    file_path = input(f"Unable to find {file_description} in current folder. Please input path of .out file for {file_description}: ")
    return file_path

def check_and_get_file_path(default_name, file_description):
    if not os.path.exists(default_name):
        return get_file_path(file_description)
    return default_name

N_out_path = check_and_get_file_path('N.out', 'N electrons state')
N_plus_1_out_path = check_and_get_file_path('N+1.out', 'N+1 electrons state')
N_minus_1_out_path = check_and_get_file_path('N-1.out', 'N-1 electrons state')

def correct_nitrogen_valence(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetExplicitValence() == 5:
            connected_atoms = [bond.GetOtherAtom(atom) for bond in atom.GetBonds()]
            oxygens = [a for a in connected_atoms if a.GetSymbol() == 'O']
            carbon = next((a for a in connected_atoms if a.GetSymbol() == 'C'), None)

            if len(oxygens) == 2 and carbon:
                bond_to_first_oxygen = mol.GetBondBetweenAtoms(atom.GetIdx(), oxygens[0].GetIdx())
                if bond_to_first_oxygen.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    bond_to_first_oxygen.SetBondType(Chem.rdchem.BondType.SINGLE)
                    oxygens[0].SetFormalCharge(-1)

                bond_to_second_oxygen = mol.GetBondBetweenAtoms(atom.GetIdx(), oxygens[1].GetIdx())
                if bond_to_second_oxygen.GetBondType() != Chem.rdchem.BondType.DOUBLE:
                    bond_to_second_oxygen.SetBondType(Chem.rdchem.BondType.DOUBLE)
                    oxygens[1].SetFormalCharge(0)

                atom.SetFormalCharge(1)

    return mol
    
def convert_to_xyz(file_path, output_xyz_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    start_index = lines.index('CARTESIAN COORDINATES (ANGSTROEM)\n') + 2

    end_index = next((i for i, line in enumerate(lines[start_index:], start=start_index) if line.strip() == ''), None)

    atom_lines = [line for line in lines[start_index:end_index] if len(line.split()) == 4]

    xyz_content = f"{len(atom_lines)}\n\n" + "".join(atom_lines)

    with open(output_xyz_path, 'w') as xyz_file:
        xyz_file.write(xyz_content)

def convert_xyz_to_mol(xyz_file_path, xyztomol_path):
    cmd = [xyztomol_path, "-ixyz", xyz_file_path, "-omol", "-m"]
    
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")

N_input = N_out_path  
xyz_file_path = 'PS.xyz'     

convert_to_xyz(N_input, xyz_file_path)

xyztomol_path = 'convert/xyztomol.exe'

convert_xyz_to_mol(xyz_file_path, xyztomol_path)
os.remove(xyz_file_path)

def find_thiophene_alpha_carbons(mol):
    thiophene_smarts = "s1cccc1"  # SMARTS pattern for thiophene
    thiophene = Chem.MolFromSmarts(thiophene_smarts)
    matches = mol.GetSubstructMatches(thiophene)

    if not matches:
        raise ValueError("No thiophene moiety found in the molecule. Script terminated.")

    alpha_carbons_pairs = []
    for match in matches:
        sulfur_index = match[0]
        sulfur_atom = mol.GetAtomWithIdx(sulfur_index)

        if sulfur_atom.GetSymbol() == 'S':
            neighbors = [neighbor.GetIdx() for neighbor in sulfur_atom.GetNeighbors() if neighbor.GetSymbol() == 'C']
            
            if len(neighbors) == 2:
                alpha_carbons_pairs.append([idx + 1 for idx in neighbors])

    return alpha_carbons_pairs

mol_file_path = 'PS.mol'
mol = Chem.MolFromMolFile(mol_file_path, sanitize=False, removeHs=False)
os.remove(mol_file_path)
mol = correct_nitrogen_valence(mol)

for atom in mol.GetAtoms():
    if atom.GetSymbol() == 'N' and atom.GetDegree() == 4:
        atom.SetFormalCharge(1)

for atom in mol.GetAtoms():
    if atom.GetSymbol() == 'N':
        double_bonds = sum(1 for bond in atom.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)      
        if double_bonds == 2:
            atom.SetFormalCharge(1)
              
for atom in mol.GetAtoms():
    if atom.GetSymbol() == 'B' and atom.GetDegree() == 4:
        atom.SetFormalCharge(-1)

for atom in mol.GetAtoms():
    if atom.GetSymbol() == 'P' and atom.GetDegree() == 6:
        atom.SetFormalCharge(-1)

for atom in mol.GetAtoms():
    if atom.GetSymbol() == 'O':
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and atom.GetDegree() == 1:
                atom.SetFormalCharge(-1)                
try:
    Chem.SanitizeMol(mol)
except Chem.rdchem.AtomValenceException as e:
    print("Error:", e)

alpha_carbons_pairs = find_thiophene_alpha_carbons(mol)

def read_file_lines(file_path: str) -> List[str]:
    try:
        with open(file_path, "r") as file:
            return file.readlines()
    except Exception as e:
        raise ValueError(f"Error reading {file_path}: {e}")

def parse_hirshfeld_data(lines: List[str], atom_indices: List[int]) -> Dict[int, float]:
    hirshfeld_data_found = False
    atom_charges = {}
    file_atom_indices = [index - 1 for index in atom_indices]

    for i in range(len(lines) - 1, -1, -1):
        line = lines[i]
        if "TOTAL" in line:
            hirshfeld_data_found = True
        elif hirshfeld_data_found and "HIRSHFELD ANALYSIS" in line:
            break
        elif hirshfeld_data_found:
            parts = line.split()
            if parts and parts[0].isdigit():
                atom_index = int(parts[0])
                if atom_index in file_atom_indices:
                    charge = float(parts[2])
                    atom_charges[atom_index] = charge
    return atom_charges

def parse_energy_value(lines: List[str]) -> float:
    for line in lines:
        if "FINAL SINGLE POINT ENERGY" in line:
            return float(line.split("FINAL SINGLE POINT ENERGY")[1].strip())
    raise ValueError("Energy value not found in file.")

def read_hirshfeld_and_energy(file_path: str, atom_indices: List[int]) -> Tuple[Dict[int, float], float]:
    lines = read_file_lines(file_path)
    atom_charges = parse_hirshfeld_data(lines, atom_indices)
    energy = parse_energy_value(lines)
    return atom_charges, energy

def process_files(alpha_carbons_pairs: List[List[int]]) -> float:
    charges_n_minus_1, energy_n_minus_1 = read_hirshfeld_and_energy(N_minus_1_out_path, [idx for pair in alpha_carbons_pairs for idx in pair])
    charges_n, energy_n = read_hirshfeld_and_energy(N_out_path, [idx for pair in alpha_carbons_pairs for idx in pair])
    charges_n_plus_1, energy_n_plus_1 = read_hirshfeld_and_energy(N_plus_1_out_path, [idx for pair in alpha_carbons_pairs for idx in pair])

    VIP = (energy_n_minus_1 - energy_n) * 27.211396641308
    VEA = (energy_n - energy_n_plus_1) * 27.211396641308
    XM = (VIP + VEA) / 2
    S = 1 / (VIP - VEA)

    Z = 0
    for pair in alpha_carbons_pairs:
        index1, index2 = pair[0] - 1, pair[1] - 1

        s_plus_1 = abs(S * (charges_n_plus_1.get(index1, 0) - charges_n.get(index1, 0)))
        s_plus_2 = abs(S * (charges_n_plus_1.get(index2, 0) - charges_n.get(index2, 0)))

        sum_s_plus = s_plus_1 + s_plus_2
        F = math.log10(abs(XM * (S + 30.75 * sum_s_plus) + (0.03822 / sum_s_plus)))
        Y = 8.612 * (charges_n.get(index1, 0) + charges_n.get(index2, 0)) - 76.0* F + 53.2494        
        Z += 10**Y
        log_k_kH = math.log10(Z)

    return log_k_kH

y_value = process_files(alpha_carbons_pairs)

def print_result(y_value):
    box_width = 41
    text = f"log(k/kH): {y_value:.2f}"
    # Calculate the spacing needed to center the text
    left_spacing = (box_width - len(text) - 2) // 2  
    right_spacing = box_width - len(text) - left_spacing - 2

    centered_line = "*" + " " * left_spacing + text + " " * right_spacing + "*"

    print("*****************************************")
    print("*                                       *")
    print("*   REACTIVITY TOWARDS SINGLET OXYGEN   *")
    print("*                                       *")
    print(centered_line)
    print("*                                       *")
    print("*****************************************")
    print("Note: This estimate is derived from a theoretical context and does not consider solvent effects.")
    print("      The reactivity is relative to the rate constant of unmodified thiophene (kH).\n")
    
print_result(y_value)



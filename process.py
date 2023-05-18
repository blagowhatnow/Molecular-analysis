#Reads smiles from a text file, assings them names, and writes the output as an sdf file

from rdkit import Chem

# Read SMILES from the text file
with open('short_listed.smi', 'r') as file:
    smiles_list = file.read().splitlines()

# Create RDKit molecules from the SMILES
molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

# Assign molecule names
for i, mol in enumerate(molecules):
    mol.SetProp("_Name", f"molecule {i+1}")

# Write molecules to an SDF file
writer = Chem.SDWriter('asinex.sdf')
for mol in molecules:
    writer.write(mol)
writer.close()








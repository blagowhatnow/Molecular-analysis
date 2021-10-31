#needs the pubchemlibrary. Use pip install pubchem 

import pubchempy

text_file = open("dataset_cleansed.smi", "r")
lines = text_file.read().splitlines()
smiles=lines

iupac_names=[]

for i in smiles:
  compounds = pubchempy.get_compounds(smiles, namespace='smiles')
  match = compounds[0]
  iupac_names.append(match.iupac_name)



print(iupac_names)

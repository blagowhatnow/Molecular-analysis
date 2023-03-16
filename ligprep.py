####DESALTS, GENERATES TAUTOMERS, PREDICTS THE IONIZATION STATE AT pH W/ DIMORPHITE DL(NEEDS TO BE DOWNLOADED) AND GENERATES CONFORMERS/ROTAMERS(A TYPE OF STEREOISOMERS)###############################


#Generate Tautomers and remove salt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolStandardize
from rdkit.Chem.SaltRemover import SaltRemover

remover = SaltRemover()

# Open the input file
with open('smiles.smi', 'r') as infile:
    # Open the output file
    with open('tautomer_smiles.smi', 'w') as outfile:
        # Loop through each line in the input file
        for line in infile:
            # Split the line into the name and SMILES string
            smiles, name = line.strip().split('\t')
            
            mol = Chem.MolFromSmiles(smiles)
            stripped_mol= remover.StripMol( mol ) 
            smiles_desalted=Chem.MolToSmiles(stripped_mol) 

            # Generate tautomers for the molecule 
            tautomers = MolStandardize.enumerate_tautomers_smiles(smiles_desalted) # Change the number of tautomers to generate
            # Write the original molecule and its tautomers to the output file
            outfile.write(name + '\t' + smiles_desalted + '\n')
            for i in range(0,len(tautomers)):
                outfile.write(name + '_tautomer_'+str(i+1)+'\t' + list(tautomers)[i] +  '\n')

import dimorphite_dl

#set these parameters for pH

min=-3.0
max=-2.0

# Open the input file
with open('tautomer_smiles.smi', 'r') as infile:
    # Open the output file
    with open('tautomer+ionized_states.smi', 'w') as outfile:
        # Loop through each line in the input file
        for line in infile:
            # Split the line into the name and SMILES string
            name, smiles = line.strip().split('\t')
            
            mols=[Chem.MolFromSmiles(smiles)]
            protonated_mols = dimorphite_dl.run_with_mol_list(mols,min_ph=5.0,max_ph=9.0,)

            outfile.write(name + '\t' + smiles + '\n')

            for i in range(0,len(protonated_mols)):
                outfile.write(name + '_ionization_state_' + str(i+1) + '\t' + Chem.MolToSmiles(protonated_mols[i]) + '\n')


# Open the input file for saving to sdf
with open('tautomer+ionized_states.smi', 'r') as infile:
    # Create a list to store the RDKit molecules
    mols = []
    # Loop through each line in the input file
    for line in infile:
        # Split the line into the SMILES string and molecule name
        name, smiles = line.strip().split('\t')
        # Create an RDKit molecule from the SMILES string
        mol = Chem.MolFromSmiles(smiles)
        # Set the molecule name
        mol.SetProp('_Name', name)
        # Generate 3D coordinates for the molecule
        AllChem.EmbedMolecule(mol)
        # Add the molecule to the list
        mols.append(mol)

# Create a multi-molecule SD writer
writer = Chem.SDWriter('output.sdf')

# Write the molecules to the SD file
for mol in mols:
    writer.write(mol)

# Close the SD writer
writer.close()

import sys
sys.path.append('/usr/local/lib/python3.7/site-packages/')

import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import subprocess

# File locations
sdfFilePath = 'output.sdf' # The input file of structures to generate conformations from
ConfoutputFilePath = 'conf_output.sdf' # Output file containing conformations for docking

inputMols = [x for x in Chem.SDMolSupplier(sdfFilePath,removeHs=False)]
# Assign atomic chirality based on the structures:
len(inputMols) # Check how many strucures


#Check that all molecules have a name
for i, mol in enumerate(inputMols):
    if mol is None:
        print('Warning: Failed to read molecule %s in %s' % (i, sdfFilePath))
    if not mol.GetProp('_Name'):

        print('Warning: No name for molecule %s in %s' % (i, sdfFilePath))

import multiprocessing

# Download this from http://pypi.python.org/pypi/futures
from concurrent import futures

# conda install progressbar
import progressbar

#Find number cores available, leave two or system might be unresponsive
numcores = multiprocessing.cpu_count()
max_workers = numcores - 2

#Knowledge based torsion generator http://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00654
# This function is called in the subprocess.
# The parameters (molecule and number of conformers) are passed via a Python

ps = AllChem.ETKDG()
ps.pruneRmsThresh=0.5
ps.numThreads=0
#Edit for number of confs desired eg n = 5

n=5

def generateconformations(m, n, name):
    m = Chem.AddHs(m)
    
    ids=AllChem.EmbedMultipleConfs(m, n, ps)
    for id in ids:
    # perform  minimization 
          AllChem.UFFOptimizeMolecule(m, confId=id)
    # EmbedMultipleConfs returns a Boost-wrapped type which)
    # cannot be pickled. Convert it to a Python list, which can.
    return m, list(ids), name


writer = Chem.SDWriter(ConfoutputFilePath)


with futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
    # Submit a set of asynchronous jobs
    jobs = []
    for mol in inputMols:
        if mol:
            name = mol.GetProp('_Name')
            job = executor.submit(generateconformations, mol, n, name)
            jobs.append(job)

    widgets = ["Generating conformations; ", progressbar.Percentage(), " ",
               progressbar.ETA(), " ", progressbar.Bar()]
    pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(jobs))
    for job in pbar(futures.as_completed(jobs)):
        mol, ids, name = job.result()
        mol.SetProp('_Name', name)
        for id in ids:
            writer.write(mol, confId=id)
writer.close()

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolStandardize

ps = AllChem.ETKDG()
ps.pruneRmsThresh = 0.5
ps.numThreads = 0
n = 5

def generateconformations(smi, n, name):
    tautomers = MolStandardize.enumerate_tautomers_smiles(smi)
    conformers = []
    for tautomer in tautomers:
        mol = Chem.AddHs(Chem.MolFromSmiles(tautomer))
        ids = AllChem.EmbedMultipleConfs(mol, n, ps)
        for id in ids:
            AllChem.UFFOptimizeMolecule(mol, confId=id)
        conformers.append((mol, list(ids)))
    return conformers, name

smi = 'CC1(C2C1C(N(C2)C(=O)C(C(C)(C)C)NC(=O)C(F)(F)F)C(=O)NC(CC3CCNC3=O)C#N)C'
sdf_output_file = 'conformers.sdf'

writer = Chem.SDWriter(sdf_output_file)

conformers, mol_name = generateconformations(smi, n, 'Molecule')
for mol, ids in conformers:
    for id in ids:
        mol_copy = Chem.Mol(mol)
        mol_copy.RemoveAllConformers()
        mol_copy.AddConformer(mol.GetConformer(id), assignId=True)
        mol_copy.SetProp('_Name', f'{mol_name}_conf{id}')
        writer.write(mol_copy)

writer.close()

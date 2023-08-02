#Kusner et al., 2017
#https://bioinformatics.stackexchange.com/questions/15645/why-do-molecular-generation-models-maximize-penalized-logp-as-a-measure-of-dru

from rdkit import Chem
from rdkit.Chem import RDConfig
import os
import sys
from rdkit.Chem import Descriptors
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

def sascore(m):
  s = sascorer.calculateScore(m)
  return s


def cycle(m,min_ring_size):
  # Define the minimum ring size you consider as "long"
  ring_info = m.GetRingInfo()
  # Count the number of rings with size greater than or equal to min_ring_size
  num_long_cycles = sum(1 for ring in ring_info.AtomRings() if len(ring) > min_ring_size)
  return num_long_cycles

def penalized_logp(m):
   logp = Descriptors.MolLogP(m)
   min_ring_size=6
   p=logp-sascore(m)-cycle(m,min_ring_size)
   return p

#Testing on Tolvaptan(Compound name)
smiles='CC1=CC=CC=C1C(=O)NC2=CC(=C(C=C2)C(=O)N3CCCC(C4=C3C=CC(=C4)Cl)O)C'
m=Chem.MolFromSmiles(smiles)
plogp=penalized_logp(m)
print(plogp)

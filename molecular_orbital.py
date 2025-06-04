import numpy as np
from sympy.physics.hydrogen import R_nl
from sympy import symbols
from sympy.utilities.lambdify import lambdify
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def show_atom_number(mol, label):
    for atom in mol.GetAtoms():
        atom.SetProp(label, str(atom.GetIdx()+1))
    return mol

target_mol_SMILE = input("")
target_mol = Chem.MolFromSmiles(target_mol_SMILE)
target_mol = Chem.AddHs(target_mol)

AllChem.EmbedMolecule(target_mol)
AllChem.MMFFOptimizeMolecule(target_mol)

conf = target_mol.GetConformer()
pos = conf.GetPositions()
show_atom_number(target_mol, 'atomNote')
Draw.MolToFile(target_mol,'testGraph.png',size=(1800,600))
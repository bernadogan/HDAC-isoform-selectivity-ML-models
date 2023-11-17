#!/usr/bin/env python

import pandas as pd
import numpy as np
from rdkit.Chem import PandasTools
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
import tmap as tm
from map4 import MAP4Calculator

"""
    Encode a molecule from a SMILES string into a fingerprint.

    Parameters
    ----------
    smiles : str
        The SMILES string defining the molecule.

    n_bits : int
        The length of the fingerprint.
    
    Returns
    -------
    array
        The fingerprint array.

    """
def smiles_to_maccs(smiles):
    # convert smiles to RDKit mol object
	mol = Chem.MolFromSmiles(smiles)
	return np.array(MACCSkeys.GenMACCSKeys(mol))
	

def smiles_to_morgan3(smiles, n_bits = 1024):
	mol = Chem.MolFromSmiles(smiles)
	return np.array(GetMorganFingerprintAsBitVect(mol, 3, nBits=n_bits))


def smiles_to_morgan3_2048(smiles, n_bits = 2048):
	mol = Chem.MolFromSmiles(smiles)
	return np.array(GetMorganFingerprintAsBitVect(mol, 3, nBits=n_bits))

def smiles_to_morganF(smiles, n_bits = 1024):
	mol = Chem.MolFromSmiles(smiles)
	return np.array(GetMorganFingerprintAsBitVect(mol, 3, nBits=n_bits, useFeatures=True))

def smiles_to_morganF_2048(smiles, n_bits = 2048):
	mol = Chem.MolFromSmiles(smiles)
	return np.array(GetMorganFingerprintAsBitVect(mol, 3, nBits=n_bits, useFeatures=True))

dim = 1024

MAP4 = MAP4Calculator(dimensions=dim)
def smiles_to_MAP4(smiles):  
    mol = Chem.MolFromSmiles(smiles)
    return np.array(MAP4.calculate(mol))

dim = 2048
MAP4_2 = MAP4Calculator(dimensions=dim)
def smiles_to_MAP4_2048(smiles):  
    mol = Chem.MolFromSmiles(smiles)
    return np.array(MAP4_2.calculate(mol))
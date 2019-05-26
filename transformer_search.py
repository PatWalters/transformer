#!/usr/bin/env python

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor, rdmolfiles, Get3DDistanceMatrix
import pandas as pd
import numpy as np


def remove_salts(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    largest = -1
    largest_frag = None
    for frag in frags:
        num_atoms = frag.GetNumAtoms()
        if num_atoms > largest:
            largest_frag = frag
    return largest_frag


def build_smiles_dictionary(df):
    smiles_dict = {}
    for (mol, name, val) in df[["Mol", "Name", "Value"]].values:
        smiles_dict[Chem.MolToSmiles(mol)] = [name, val]
    return smiles_dict


def clear_atom_maps(mol):
    for atm in mol.GetAtoms():
        atm.SetAtomMapNum(0)


def get_bond_lengths(mol):
    """
    Return a list of bond lengths
    :param mol: input molecule
    :return: list of bond lengths
    """
    dm = Get3DDistanceMatrix(mol)
    bnd_list = []
    for bnd in mol.GetBonds():
        start = bnd.GetBeginAtomIdx()
        end = bnd.GetEndAtomIdx()
        bnd_list.append(dm[start, end])
    return bnd_list


def scale_molecule(mol, factor=1.5):
    """
    Scale the bond lengths in a molecule
    :param mol: input molecule
    :param factor: scaling factor
    :return: None
    """
    mean_dist = np.mean(get_bond_lengths(mol))
    factor = factor / mean_dist
    matrix = np.zeros((4, 4), np.float)
    for i in range(3):
        matrix[i, i] = factor
        matrix[3, 3] = 1
    AllChem.TransformMol(mol, matrix)


def run_transforms(rxn_file, data_file):
    rxn = AllChem.ReactionFromRxnFile(rxn_file)
    reactant_template = Chem.Mol(rxn.GetReactantTemplate(0))
    product_template = Chem.Mol(rxn.GetProductTemplate(0))
    scale_molecule(reactant_template)
    scale_molecule(product_template)
    clear_atom_maps(reactant_template)
    clear_atom_maps(product_template)
    df = pd.read_csv(data_file)
    mol_list = [Chem.MolFromSmiles(x) for x in df.SMILES]
    df["Mol"] = [remove_salts(x) for x in mol_list]
    smiles_dict = build_smiles_dictionary(df)
    output_list = []
    used = set()
    for (mol, name, val) in df[["Mol", "Name", "Value"]].values:
        prods = rxn.RunReactants([mol])
        for prod in prods:
            prod_mol = prod[0]
            prod_smiles = Chem.MolToSmiles(prod_mol)
            prod_lookup = smiles_dict.get(prod_smiles)
            if prod_lookup is not None:
                if prod_smiles not in used:
                    rdDepictor.GenerateDepictionMatching2DStructure(mol, reactant_template)
                    output_list.append([mol, name, val])
                    prod_name, prod_val = prod_lookup
                    prod_mol.UpdatePropertyCache()
                    rdDepictor.GenerateDepictionMatching2DStructure(prod_mol, product_template)
                    output_list.append([prod_mol, prod_name, prod_val])
                    used.add(prod_smiles)
    return output_list


def save_output(output_list, file_name):
    writer = rdmolfiles.SDWriter(file_name)
    for mol, name, val in output_list:
        mol.SetProp("_Name", name)
        mol.SetProp("Value", str(val))
        writer.write(mol)


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} rxn.rxn infile.csv outfile.sdf")
        sys.exit(1)
    out_list = run_transforms(sys.argv[1], sys.argv[2])
    save_output(out_list, sys.argv[3])


main()

#!/usr/bin/env python

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor, rdmolfiles, Get3DDistanceMatrix
import pandas as pd
import numpy as np
from docopt import docopt

doc = """Usage:
transformer_search.py --rxn REACTION_FILE --in CSV_FILE --out OUTPUT_FILE[--smarts SMARTS_PATTERN]

Options: 
--rxn REACTION_FILE reaction transform as an MDL rxn file
--in CSV_FILE  input data, columns are SMILES, Name, Value
--out OUTPUT_FILE sdf file for output data
--smarts SMARTS_PATTERN optional smarts pattern that all reactants must match
"""


def remove_salts(mol):
    """
    Remove everything but the largest fragment
    :param mol: input molecule
    :return: largest fragment as a molecule
    """
    frags = Chem.GetMolFrags(mol, asMols=True)
    largest = -1
    largest_frag = None
    for frag in frags:
        num_atoms = frag.GetNumAtoms()
        if num_atoms > largest:
            largest_frag = frag
            largest = num_atoms
    return largest_frag


def build_smiles_dictionary(df):
    """
    Convert input dataframe to a dictionary
    :param df: input dataframe
    :return: output dictionary
    """
    smiles_dict = {}
    for (mol, name, val) in df[["Mol", "Name", "Value"]].values:
        smiles_dict[Chem.MolToSmiles(mol)] = [name, val]
    return smiles_dict


def clear_atom_maps(mol):
    """
    Set all atom maps to 0
    :param mol:  input molecule
    :return: None
    """
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


def run_transforms(rxn_file, data_file, smarts):
    """
    Given a reaction and a file with SMILES and activity values, generate products and check whether products are in
    the input
    :param rxn_file: reaction file name
    :param data_file: data file name
    :param smarts: smarts that must be matched by input molecules
    :return: list of output molecules
    """
    if smarts:
        smarts_mol = Chem.MolFromSmarts(smarts)
        if smarts_mol is None:
            print(f"Could not parse SMARTS {smarts}")
            sys.exit(1)
    rxn = AllChem.ReactionFromRxnFile(rxn_file)
    reactant_template = Chem.Mol(rxn.GetReactantTemplate(0))
    product_template = Chem.Mol(rxn.GetProductTemplate(0))
    # In order to make the output easier to understand, we will align the molecules to the reaction
    # The aligned structures will look funny unless we scale the reaction and product templates
    scale_molecule(reactant_template)
    scale_molecule(product_template)
    # Clear the atom maps from the templates, may not be necessary, but can't hurt
    clear_atom_maps(reactant_template)
    clear_atom_maps(product_template)
    # Read in the data
    df = pd.read_csv(data_file)
    # Create a list of input molecules
    mol_list = [Chem.MolFromSmiles(x) for x in df.SMILES]
    df["Mol"] = [remove_salts(x) for x in mol_list]
    # Build a dictionary of the input SMILES
    smiles_dict = build_smiles_dictionary(df)
    output_list = []
    used = set()
    # Loop over the molecules, apply the reaction and check if the product SMILES is in smiles_dict
    # If it is, we have a pair
    for (mol, name, val) in df[["Mol", "Name", "Value"]].values:
        if smarts and not mol.HasSubstructMatch(smarts_mol):
            continue
        prods = rxn.RunReactants([mol])
        for prod in prods:
            prod_mol = prod[0]
            prod_smiles = Chem.MolToSmiles(prod_mol)
            prod_lookup = smiles_dict.get(prod_smiles)
            if prod_lookup is not None:
                # skip duplicate products
                if prod_smiles not in used:
                    # Generate a depiction aligned to the reactant template
                    rdDepictor.GenerateDepictionMatching2DStructure(mol, reactant_template)
                    output_list.append([mol, name, val])
                    prod_name, prod_val = prod_lookup
                    prod_mol.UpdatePropertyCache()
                    # Generate a depiction aligned to the product template
                    rdDepictor.GenerateDepictionMatching2DStructure(prod_mol, product_template)
                    output_list.append([prod_mol, prod_name, prod_val])
                    used.add(prod_smiles)
    return output_list


def save_output(output_list, file_name):
    """
    Write the output to an SD file
    :param output_list: list of pairs
    :param file_name: output file name
    :return: None
    """
    writer = rdmolfiles.SDWriter(file_name)
    for mol, name, val in output_list:
        mol.SetProp("_Name", name)
        mol.SetProp("Value", str(val))
        writer.write(mol)


def main():
    cmd_input = docopt(doc)
    rxn_file_name = cmd_input.get("--rxn")
    csv_file_name = cmd_input.get("--in")
    out_file_name = cmd_input.get("--out")
    smarts = cmd_input.get("--smarts")
    out_list = run_transforms(rxn_file_name, csv_file_name, smarts)
    save_output(out_list, out_file_name)


main()

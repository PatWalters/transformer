#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor, TemplateAlign, rdmolfiles
import pandas as pd
import sys


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


def run_transforms(rxn_file, data_file):
    rxn = AllChem.ReactionFromRxnFile(rxn_file)
    reactant_template = Chem.Mol(rxn.GetReactantTemplate(0))
    product_template = Chem.Mol(rxn.GetProductTemplate(0))
    rdDepictor.Compute2DCoords(reactant_template)
    rdDepictor.Compute2DCoords(product_template)
    clear_atom_maps(reactant_template)
    clear_atom_maps(product_template)
    df = pd.read_csv(data_file)
    mol_list = [Chem.MolFromSmiles(x) for x in df.SMILES]
    df["Mol"] = [remove_salts(x) for x in mol_list]
    smiles_dict = build_smiles_dictionary(df)
    output_list = []
    for (mol, name, val) in df[["Mol", "Name", "Value"]].values:
        prods = rxn.RunReactants([mol])
        for prod in prods:
            prod_mol = prod[0]
            prod_smiles = Chem.MolToSmiles(prod_mol)
            prod_lookup = smiles_dict.get(prod_smiles)
            if prod_lookup is not None:
                print(Chem.MolToSmiles(prod_mol))
                rdDepictor.GenerateDepictionMatching2DStructure(mol,reactant_template)
                output_list.append([mol, name, val])

                prod_name, prod_val = prod_lookup
                rdDepictor.GenerateDepictionMatching2DStructure(prod_mol, product_template)
                output_list.append([prod_mol, prod_name, prod_val])
    return output_list


def save_output(output_list, file_name):
    writer = rdmolfiles.SDWriter(file_name)
    for mol, name, val in output_list:
        mol.SetProp("_Name",name)
        mol.SetProp("Value",str(val))
        writer.write(mol)


def main():
    out_list = run_transforms("test2.rxn", "CHEMBL1949661.csv")
    save_output(out_list,"out.sdf")


def test():
    mol = Chem.MolFromSmiles("FC(F)(F)Oc1cccc(-n2nnc3ccc(NC4CCOCC4)nc32)c1")
    rdDepictor.Compute2DCoords(mol)
    tmplt = Chem.MolFromSmarts("c1cnc2nnnc2c1")
    rdDepictor.Compute2DCoords(tmplt)
    TemplateAlign.AlignMolToTemplate2D(mol, tmplt, clearConfs=True)
    print(mol.HasSubstructMatch(tmplt))

main()



#!/usr/bin/env python

import sys
from rdkit import Chem
from rdkit.Chem import AllChem

# simple reaction script for debugging

if len(sys.argv) != 3:
    print(f"usage: {sys.argv[0]} rxn.rxn infile.smi")
    sys.exit(0)

rxn_file_name = sys.argv[1]
rxn = AllChem.ReactionFromRxnFile(rxn_file_name)
if rxn is None:
    print(f"Reaction {rxn_file_name} is not valid")
    sys.exit(1)
suppl = Chem.SmilesMolSupplier(sys.argv[2])

for mol in suppl:
    if mol is None:
        continue
    prods = rxn.RunReactants([mol])
    if len(prods):
        input_smiles = Chem.MolToSmiles(mol)
        name = mol.GetProp("_Name")
        smiles_list = [Chem.MolToSmiles(x[0]) for x in prods]
        smiles_list = list(set(smiles_list))
        for smiles in smiles_list:
            print(input_smiles, name)
            print(smiles, name)

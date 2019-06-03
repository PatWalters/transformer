This code accompanys my blog post [Using Reaction Transforms to Understand SAR](https://practicalcheminformatics.blogspot.com/2019/06/using-reaction-transforms-to-understand.html) 

```
Usage:
transformer_search.py --rxn REACTION_FILE --in CSV_FILE --out OUTPUT_FILE[--smarts SMARTS_PATTERN]

Options:
--rxn REACTION_FILE reaction transform as an MDL rxn file
--in CSV_FILE  input data, columns are SMILES, Name, Value
--out OUTPUT_FILE sdf file for output data
--smarts SMARTS_PATTERN optional smarts pattern that all reactants must match
```
The "--smarts" argument won't typically be required.  The only time you may want to use this is when you want to exclude certain molecules.  For instance, let's say that you want molecules with methyl attached to the scaffold, but not molecules with ethyl.  If you use this argument, only molecules matching the SMARTS will be processed by the reaction transform. 


### Installation
This script requires [the RDKit](https://www.rdkit.org/docs/Install.html).  The included Jupyter notebook also requires the [Seaborn graphics library](https://seaborn.pydata.org/).
```
pip install seaborn
```
### Example
To run the included example:
```
transformer_search.py --rxn rxn.rxn --in CHEMBL1949661.csv --out out.sdf
```

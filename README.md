Using reaction transforms to identify spectific matched molecular pairs (MMPs).  This code accompanys my blog post [Using Reaction Transforms to Understand SAR](https://practicalcheminformatics.blogspot.com/2019/06/using-reaction-transforms-to-understand.html) 

```
Usage:
transformer_search.py --rxn REACTION_FILE --in CSV_FILE --out OUTPUT_FILE[--smarts SMARTS_PATTERN]

Options:
--rxn REACTION_FILE reaction transform as an MDL rxn file
--in CSV_FILE  input data, columns are SMILES, Name, Value
--out OUTPUT_FILE sdf file for output data
--smarts SMARTS_PATTERN optional smarts pattern that all reactants must match
```


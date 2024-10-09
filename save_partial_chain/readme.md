Extract part of the protein chain from a PDB file (Informatics)

This tool is used for extracting certain residues in a protein structure, and save a new PDB file with the extracted residues. The tool is helpful for extracting a part of a protein chain, for example, a domain, a binding site, or a loop region.

# Input
input_pdb: The input PDB file to extract substructure from.
retained_part: A dictionary with the residues to keep in each chain. The dictionary keys are the chain IDs, and the values are strings. A range of residues can be represented with columns, i.e. '5:20' will keep residues 5 to 20 (inclusive for both ends). Additional residues can be added with comma, like '10,15,18'. As an example, {'A': '1:5', 'B': '3:6,9,11'} will keep residues 1, 2, 3, 4 and 5 in chains A and residues 3, 4, 5, 6, 9 and 11 in B.

# Output
A pdb with only the specified residues retained.


## Install
### requirements
biopython
pdb-tools
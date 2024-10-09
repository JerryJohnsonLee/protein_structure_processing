Cap N and C terminal of a protein chain (Modelling)

Add ACE cap to the N-terminal and NME cap to the C-terminal of a protein chain. The input structure should contain only one chain. When the input structure contains hydrogens, the added caps will also contain hydrogens. Otherwise the caps will not have hydrogens.

# Input
input_pdb - file path to the single-chain PDB structure that needs to be capped
N_terminal_cap - either "ACE" or "none". If "none", then the N terminal will not be capped
C_terminal_cap - either "NME" or "none". If "none", then the C terminal will not be capped

# Output
A PDB file with the appropriate terminal caps added

## Install
### requirements
biopython
numpy
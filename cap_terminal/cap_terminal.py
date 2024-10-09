from Bio import PDB
import os
from typing import Literal, Tuple
import numpy as np

parser = PDB.PDBParser(QUIET=True)

def standard_svd(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Perform singular value decomposition on the input matrix. 
    # The singular values are provided as a diagonal matrix.
    u, s, v = np.linalg.svd(matrix)
    S = np.zeros((u.shape[1],v.shape[0]))
    S[:len(s), :len(s)] = np.diag(s)
    return u, S, v

def find_rotation(a: np.ndarray, b: np.ndarray):
    # Compute RMSD using the Kabsch algorithm

    dot_prod = np.dot(a.T, b)
    v, s, wt = standard_svd(dot_prod)
    d = np.sign(np.linalg.det(np.dot(v, wt).T))
    e = np.eye(3)
    e[2,2] = d
    u = np.dot(np.dot(wt.T, e), v.T)     # rotation matrix
    return u

def find_H_position(N_pos, CA_pos, C_pos):
    d_N_CA = np.linalg.norm(CA_pos - N_pos)
    d_CA_C = np.linalg.norm(C_pos - CA_pos)
    angle_N_CA_C = np.arccos(np.dot(N_pos - CA_pos, C_pos - CA_pos) / (d_N_CA * d_CA_C))
    cos_supp = np.cos(np.pi - angle_N_CA_C)
    tan_supp = np.tan(np.pi - angle_N_CA_C)
    reference_pos = CA_pos + (CA_pos - C_pos) / d_CA_C * d_N_CA * np.sqrt(3) / (cos_supp * tan_supp * (np.sqrt(3) / tan_supp + 1))
    H_pos = N_pos + (N_pos - reference_pos) / np.linalg.norm(reference_pos - N_pos) * 1.01
    return H_pos

def run(input_pdb, N_terminal_cap: Literal["ACE", "CH3", "none"], C_terminal_cap: Literal["NH2", "NME", "CH3", "OH", "none"]):
    '''
    use OpemMM pdbfixer to add terminal caps to a protein structure
    '''
    
    # Load the PDB structure
    structure = parser.get_structure(input_pdb, input_pdb)

    # make sure there is only one chain, and there is no gap in the chain
    assert len(structure) == 1 and len(list(structure[0].get_chains())) == 1, \
        "The input PDB file should contain 1 model and 1 chain"
    chain = list(structure[0].get_chains())[0]
    residue_ids = np.array([residue.id[1] for residue in chain])
    # assert (residue_ids[1:] - residue_ids[:-1] == 1).all(), "The chain should not have any gaps"

    # check whether the original structure contains hydrogens
    has_hydrogens = any([atom.element == 'H' for atom in chain.get_atoms()])

    # process N terminal
    N_term_id = residue_ids[0]
    N_term_residue = chain[int(N_term_id)]

    # make sure required atoms are present in the N terminal residue
    assert 'N' in N_term_residue, "Missing N atom in the N terminal residue"
    assert 'CA' in N_term_residue, "Missing CA atom in the N terminal residue"
    if has_hydrogens:
        assert 'H' in N_term_residue, "Missing H atom in the N terminal residue"
    else:
        assert 'C' in N_term_residue, "Missing C atom in the N terminal residue"

    N_term_N_pos = N_term_residue['N'].coord
    N_term_CA_pos = N_term_residue['CA'].coord
    if has_hydrogens:
        N_term_H_pos = N_term_residue['H'].coord
    else:
        N_term_C_pos = N_term_residue['C'].coord
        N_term_H_pos = find_H_position(N_term_N_pos, N_term_CA_pos, N_term_C_pos)

    # read in template structure for N terminal cap
    if N_terminal_cap != "none":
        if not os.path.exists(f"cap_template/{N_terminal_cap}.pdb"):
            raise FileNotFoundError(f"Template structure for {N_terminal_cap} not found")
        N_cap_tmeplate = parser.get_structure(N_terminal_cap, f"cap_template/{N_terminal_cap}.pdb")
        reference_coords = np.array([N_term_H_pos, N_term_CA_pos]) - N_term_N_pos
        template_coords = np.array([atom.coord for atom in N_cap_tmeplate.get_atoms()])
        rotation = find_rotation(np.array([[0, -1, 0], [0.866, 0.5, 0]]), reference_coords)
        new_coords = np.dot(template_coords, rotation.T) + N_term_N_pos
        
        # prepare cap as a new residue
        N_term_cap_residue = PDB.Residue.Residue((' ', N_term_id - 1, ' '), N_terminal_cap, ' ')
        for atm, coord in zip(N_cap_tmeplate.get_atoms(), new_coords):
            if atm.name in ["N", "H", "H1"]:
                continue
            if not has_hydrogens and atm.element == "H":
                continue
            new_atom = PDB.Atom.Atom(atm.name, coord, atm.bfactor, atm.occupancy, atm.altloc, atm.fullname, atm.serial_number, atm.element)
            N_term_cap_residue.add(new_atom)
        
        # insert the cap residue to the beginning of the chain
        chain.insert(0, N_term_cap_residue)


    # process C terminal
    C_term_id = residue_ids[-1]
    C_term_residue = chain[int(C_term_id)]

    # make sure required atoms are present in the C terminal residue
    assert 'C' in C_term_residue, "Missing C atom in the C terminal residue"
    assert 'CA' in C_term_residue, "Missing CA atom in the C terminal residue"
    assert 'O' in C_term_residue, "Missing O atom in the C terminal residue"

    C_term_C_pos = C_term_residue['C'].coord
    C_term_CA_pos = C_term_residue['CA'].coord
    C_term_O_pos = C_term_residue['O'].coord

    # read in template structure for C terminal cap
    if C_terminal_cap != "none":
        if not os.path.exists(f"cap_template/{C_terminal_cap}.pdb"):
            raise FileNotFoundError(f"Template structure for {C_terminal_cap} not found")
        C_cap_template = parser.get_structure(C_terminal_cap, f"cap_template/{C_terminal_cap}.pdb")
        reference_coords = np.array([C_term_O_pos, C_term_CA_pos]) - C_term_C_pos
        template_coords = np.array([atom.coord for atom in C_cap_template.get_atoms()])
        rotation = find_rotation(np.array([[0, 1.26, 0], [-0.866, -0.5, 0]]), reference_coords)
        new_coords = np.dot(template_coords, rotation.T) + C_term_C_pos

        # prepare cap as a new residue
        C_term_cap_residue = PDB.Residue.Residue((' ', C_term_id + 1, ' '), C_terminal_cap, ' ')
        for atm, coord in zip(C_cap_template.get_atoms(), new_coords):
            if atm.name in ["C", "O", "HC"]:
                continue
            if not has_hydrogens and atm.element == "H":
                continue
            new_atom = PDB.Atom.Atom(atm.name, coord, atm.bfactor, atm.occupancy, atm.altloc, atm.fullname, atm.serial_number, atm.element)
            C_term_cap_residue.add(new_atom)
        
        # add the cap residue to the end of the chain
        chain.add(C_term_cap_residue)



    # save final structure
    output_name = input_pdb.replace(".pdb", f".capped_{N_terminal_cap}_{C_terminal_cap}.pdb")
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_name)
    print("Capped structure saved to", output_name)



if __name__ == "__main__":
    run('single_chain_protonated.pdb', 'abc', 'CH3')
    # run('4DJH.pdb', 'ACE', 'NME')

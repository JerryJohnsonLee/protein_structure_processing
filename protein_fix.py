import os
import shutil
from datetime import datetime

from Bio import Align
from Bio.PDB import PDBParser
from Bio.Align import substitution_matrices

import warnings
warnings.filterwarnings("ignore", message=".*discontinuous at line.*")

PROMOD_EXECUTABLE = os.path.join(os.environ['PROMOD3_ROOT'], 'src', "promod-3.4.1", "build", "stage", "bin", "pm")
parser = PDBParser()

resname_mapping = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C', 
    'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P', 
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def create_timestamp() -> str:
  # helper function to create a unique timestamp
  dt = str(datetime.now())
  return dt.replace("-", "_").replace(":", "_").replace(" ", "_")


def get_structure_sequence(pdb_path):
    '''
    use biopython to get the sequence of the structure that has 3d coordinates
    returns a dictionary with chain id as key and sequence as value
    '''
    seq_dict = {}
    pdb_name = os.path.basename(pdb_path)

    structure = parser.get_structure(pdb_name, pdb_path)
    model = structure[0]
    for chain in model:
        standard_residues = [res for res in chain if res.get_resname() in resname_mapping]
        if len(standard_residues) == 0:
            continue
        seq = ''.join([resname_mapping[res.get_resname()] for res in standard_residues])
        seq_dict[chain.id] = seq
    return seq_dict

def get_full_sequence(pdb_path):
    '''
    get the full sequence of the structure from the header of the pdb file
    returns a dictionary with chain id as key and sequence as value
    '''
    seq_dict = {}
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('SEQRES'):
                chain_id = line[11]
                seq_3_letters = line[19:].strip().split()
                seq = ''.join([resname_mapping[res] for res in seq_3_letters])
                if chain_id in seq_dict:
                    seq_dict[chain_id] += seq
                else:
                    seq_dict[chain_id] = seq
    return seq_dict

def align_sequence(full_seq, structure_seq, ignore_terminal_gaps=True):
    # use biopython to align two sequences
    aligner = Align.PairwiseAligner()
    matrix = substitution_matrices.load("BLOSUM62")
    for i in range(len(str(matrix.alphabet))):
        res = matrix.alphabet[i]
        matrix['X'][res] = 0
        matrix[res]['X'] = 0
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    alignments = aligner.align(full_seq, structure_seq)
    alignment = str(alignments[0]).split("\n")
    aligned_full_seq = alignment[0]
    aligned_structure_seq = alignment[2]

    if ignore_terminal_gaps:
    # remove terminal gaps from the alignment
        while aligned_structure_seq[0] == '-':
            aligned_structure_seq = aligned_structure_seq[1:]
            aligned_full_seq = aligned_full_seq[1:]
        while aligned_structure_seq[-1] == '-':
            aligned_structure_seq = aligned_structure_seq[:-1]
            aligned_full_seq = aligned_full_seq[:-1]
    return aligned_full_seq, aligned_structure_seq


def align_sequence_dict(full_seq_dict, structure_seq_dict, workdir):
    '''
    align the full sequences from the header and the 3d structure
    returns a dictionary with chain id as key and a tuple of the aligned sequences as value
    '''
    alignment_files = []
    assert set(full_seq_dict.keys()) == set(structure_seq_dict.keys())
    for chain_id in full_seq_dict:
        full_seq = full_seq_dict[chain_id]
        structure_seq = structure_seq_dict[chain_id]

        aligned_full_seq, aligned_structure_seq = align_sequence(full_seq, structure_seq)
        alignment_file = os.path.join(workdir, f'{chain_id}_aligned.fasta')
        with open(alignment_file, 'w') as f:
            f.write(f'>Full sequence|Chain {chain_id}\n{aligned_full_seq}\n')
            f.write(f'>Structure sequence|Chain {chain_id}\n{aligned_structure_seq}\n')
        alignment_files.append(alignment_file)
    return alignment_files       


        
        

def run(input_pdb: str):
    """
    Checks a PDB file for common errors and fixes them.
    :param input_pdb: Path to the input PDB file.
    :return: Path to the adjusted PDB file.
    """

    structure_seq_dict = get_structure_sequence(input_pdb)
    full_seq_dict = get_full_sequence(input_pdb)
    timestamp = create_timestamp()
    workdir = timestamp
    output_file = input_pdb.replace('.pdb', f'_fixed.{timestamp}.pdb')
    os.makedirs(workdir)
    alignment_files = align_sequence_dict(full_seq_dict, structure_seq_dict, workdir)
    os.system(f'{PROMOD_EXECUTABLE} homology_modelling.py {input_pdb} {output_file} {" ".join(alignment_files)}')
    shutil.rmtree(workdir)
    print("Fixed structure saved at", output_file)


if __name__ == "__main__":
    # input_pdb = sys.argv[1]
    input_pdb = "4DJH.pdb"
    run(input_pdb)
from datetime import datetime
import os
import shutil

from Bio import PDB


def merge_pdb_files(pdb_file_list, output_file, pdb_name):
    # Create a PDB structure builder
    parser = PDB.PDBParser(QUIET=True)
    io = PDB.PDBIO()
    
    # Create an empty structure and a single model
    merged_structure = PDB.Structure.Structure(pdb_name)
    model = PDB.Model.Model(0)  # Single model with ID 0
    merged_structure.add(model)
    
    # Chain ID management (to avoid duplication)
    chain_id_map = {}
    
    for pdb_file in pdb_file_list:
        # Parse the PDB file
        structure = parser.get_structure(pdb_file, pdb_file)
        
        for chain in structure.get_chains():
            # Create a unique chain ID if necessary
            if chain.id in chain_id_map:
                chain_id_map[chain.id] += 1
                new_chain_id = f"{chain.id}_{chain_id_map[chain.id]}"
            else:
                chain_id_map[chain.id] = 0
                new_chain_id = chain.id
            
            new_chain = PDB.Chain.Chain(new_chain_id)
            
            for res in chain.get_residues():
                new_residue = PDB.Residue.Residue(res.id, res.resname, res.segid)
                new_chain.add(new_residue)
                
                for atom in res.get_atoms():
                    new_atom = PDB.Atom.Atom(
                        atom.name,
                        atom.coord,
                        atom.bfactor,
                        atom.occupancy,
                        atom.altloc,
                        atom.fullname,
                        atom.serial_number,
                        atom.element,
                    )
                    new_residue.add(new_atom)
            
            model.add(new_chain)

    # Write the merged structure to an output file
    io.set_structure(merged_structure)
    io.save(output_file)

def create_timestamp() -> str:
  # helper function to create a unique timestamp
  dt = str(datetime.now())
  return dt.replace("-", "_").replace(":", "_").replace(" ", "_")

def run(input_pdb, retained_part: dict):
    timetag = create_timestamp()
    workdir = "./" + timetag
    os.makedirs(workdir)
    shutil.copy(input_pdb, workdir)

    # use pdbtools to save all single chains
    for chain in retained_part:
        os.system(f"pdb_selchain -{chain} {workdir}/{input_pdb} > {workdir}/chain_{chain}.pdb")
        if retained_part[chain] in ["all", ""]:
            shutil.copy(f"{workdir}/chain_{chain}.pdb", f"{workdir}/chain_{chain}_selected.pdb")
            print("All residues in chain {} are retained.".format(chain))
        else:
            os.system(f"pdb_selres -{retained_part[chain]} {workdir}/chain_{chain}.pdb > {workdir}/chain_{chain}_selected.pdb")

    # merge all selected chains using biopython
    final_pdb = input_pdb.replace(".pdb", f".selected.{timetag}.pdb")
    merge_pdb_files([f"{workdir}/chain_{chain}_selected.pdb" for chain in retained_part], final_pdb, input_pdb.replace(".pdb", ""))
    print(f"Processed PDB file is: {final_pdb}")
          

if __name__ == "__main__":
    input_pdb = "4DJH.pdb"
    retained_part = {
        "A": "55:347",
    }
    run(input_pdb, retained_part)
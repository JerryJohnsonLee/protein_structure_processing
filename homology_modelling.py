# run through the promod3 docker
import os
import sys
from ost import io, seq
from promod3 import modelling, loop


def main(template_pdb, output_pdb, aligned_chains):
    '''
    use promod3 to model a multi-chain protein structure using a template structure
    structure as the reference
    '''
    # setup
    merge_distance = 4
    fragment_db = loop.LoadFragDB()
    structure_db = loop.LoadStructureDB()
    torsion_sampler = loop.LoadTorsionSamplerCoil()

    # read in template model
    template = io.LoadPDB(template_pdb)

    alignment_lst = seq.AlignmentList()
    for alignment in aligned_chains:
        chain_id = os.path.basename(alignment)[0]
        alignment_obj = io.LoadAlignment(alignment)
        alignment_obj.AttachView(1, template.Select('chain=' + chain_id))
        alignment_lst.append(alignment_obj)

    # create raw model from alignment
    mhandle = modelling.BuildRawModel(alignment_lst)

    # terminal gaps can be ignored
    modelling.RemoveTerminalGaps(mhandle)

    # perform loop modelling to close all gaps
    modelling.CloseGaps(mhandle, merge_distance, fragment_db,
                        structure_db, torsion_sampler)

    # build sidechains
    modelling.BuildSidechains(mhandle, merge_distance, fragment_db,
                            structure_db, torsion_sampler)

    # minimize energy of final model using molecular mechanics
    modelling.MinimizeModelEnergy(mhandle)

    # check final model and report issues
    modelling.CheckFinalModel(mhandle)

    # extract final model
    final_model = mhandle.model
    io.SavePDB(final_model, output_pdb)

    print("Modelling complete.")


if __name__ == '__main__':
    template_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    aligned_chains = sys.argv[3:]
    main(template_pdb, output_pdb, aligned_chains)
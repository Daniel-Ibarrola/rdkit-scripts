from plip.structure.preparation import PDBComplex


def extract_all_ligands_plip(file_name):
    """ Extract the ligand from a pdb using plip. """
    pdb = PDBComplex()
    pdb.load_pdb(file_name)
    pdb.analyze()

    ligands = {}
    for lig in pdb.ligands:
        ligand_id = ":".join([lig.hetid, lig.chain, str(lig.position)])
        ligands[ligand_id] = lig.mol

    return ligands



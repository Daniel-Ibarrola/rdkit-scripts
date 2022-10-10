import rdkit.Chem
from rdkit import Chem

# TODO: find the names of other ions and small molecules
non_ligand_residues = [
    "HOH",
    "DOD",
]


def is_ligand(residue):
    """ Returns whether a residue is a ligand.

        Parameters
        ----------
        residue : rdkit.Chem.AtomPDBResidueInfo
            The residue info

        Returns
        -------
        bool
    """
    if residue.GetIsHeteroAtom() and \
            residue.GetResidueName() not in non_ligand_residues:
        return True
    return False


def get_ligand_names(pdb_file):
    """ Get the names of the ligands on a pdb.

       Parameters
       ----------
        pdb_file : str
            Path to a pdb file.

       Returns
       -------
       list[str]
            List with the ligand ids
    """
    pl_complex: rdkit.Chem.Mol
    atom: rdkit.Chem.Atom
    residue_info: rdkit.Chem.AtomPDBResidueInfo

    pl_complex = Chem.MolFromPDBFile(pdb_file)
    ligands = {}  # Maps residues to number of atoms

    for atom in pl_complex.GetAtoms():
        residue_info = atom.GetPDBResidueInfo()
        residue_name = residue_info.GetResidueName()
        if is_ligand(residue_info):
            try:
                ligands[residue_name] += 1
            except KeyError:
                ligands[residue_name] = 1

    # Returns molecules with more than 5 atoms
    return [lig for lig in ligands.keys() if ligands[lig] > 5]


def get_ligand(pl_complex, ligand_id):
    """ Extract ligand from protein-ligand complex.

        Parameters
        ----------
        pl_complex : rdkit.Mol

        ligand_id : str

        Returns
        -------
        ligand : rdkit.RWMol
    """
    # Find ligand
    ligand = Chem.RWMol()
    atom_index_mapper = {}
    atom_bonds = []

    # First add the atoms
    for atom in pl_complex.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info.GetResidueName() == ligand_id:
            new_index = ligand.AddAtom(atom)
            atom_index_mapper[atom.GetIdx()] = new_index
            atom_bonds.append(atom.GetBonds())

    # TODO: add bond metadata
    # Now add the bonds
    for bonds_ in atom_bonds:
        for bond in bonds_:
            begin_atom = atom_index_mapper[bond.GetBeginAtomIdx()]
            end_atom = atom_index_mapper[bond.GetEndAtomIdx()]
            # Avoid repeated bonds. The graph is directed so if end atom > begin atom
            # the bond has already been processed
            if end_atom < begin_atom:
                ligand.AddBond(begin_atom, end_atom)

    return ligand


def split_protein_ligand_complex(pdb_file):
    """ Takes a pdb file and splits it into protein and ligand.

        Parameters
        ----------
        pdb_file : str
            Path to a pdb file.

        Returns
        -------
        protein : rdkit.RWMol
            The protein.

        ligand : rdkit.RWMol
            The ligand.
    """
    pl_complex: rdkit.Chem.Mol
    atom: rdkit.Chem.Atom
    residue_info: rdkit.Chem.AtomPDBResidueInfo

    pl_complex = Chem.MolFromPDBFile(pdb_file)


if __name__ == "__main__":
    test_file_1 = "./test_cases/eralpha/1qku/1qku.pdb"
    test_file_2 = "./test_cases/rhinovirus/1ncr/1ncr.pdb"

    pl_complex_1 = Chem.MolFromPDBFile(test_file_1)
    ligand_1 = get_ligand(pl_complex_1, "EST")
    print(ligand_1.GetNumAtoms())

    pl_complex_2 = Chem.MolFromPDBFile(test_file_2)
    ligand_2 = get_ligand(pl_complex_2, "W11")
    print(ligand_2.GetNumAtoms())

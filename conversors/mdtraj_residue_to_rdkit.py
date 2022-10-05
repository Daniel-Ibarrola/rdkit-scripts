import mdtraj as mdt
from rdkit import Chem


def residue_to_rdkit(residue_index, topology):
    """ Transform a residue into a rdkit molecule.

        Parameters
        ----------
        residue_index : int
            The index of the residue in the topology

        topology : mdt.Topology
            A topology object

        Returns
        -------
        molecule : rdkit.RWMol
            An editable molecule
    """

    molecule = Chem.RWMol()
    residue_atoms = []
    atom_index_mapper = {}

    for atom in topology.residue(residue_index).atoms:
        index_topology = atom.index
        element = atom.element.symbol
        if element == "D":
            element = "H"

        index_rd_mol = molecule.AddAtom(Chem.Atom(element))
        residue_atoms.append(index_topology)
        atom_index_mapper[index_topology] = index_rd_mol

    bonds = []
    for atom_1, atom_2 in topology.bonds:
        if atom_1.index in residue_atoms and atom_2.index in residue_atoms:
            index_1 = atom_index_mapper[atom_1.index]
            index_2 = atom_index_mapper[atom_2.index]
            if (index_1, index_2) not in bonds:
                molecule.AddBond(index_1, index_2)
                bonds.append((index_1, index_2))

    return molecule

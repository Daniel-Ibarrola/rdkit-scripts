from rdkit import Chem


def topology_to_rdkit(topology):
    """ Create and rdkit molecule from a topology.

        Parameters
        ----------
        topology : mdtraj.Topology
            A topology object

        Returns
        -------
        rdkit.RWMol
            An editable molecule.
    """
    molecule = Chem.RWMol()
    for atom in topology.atoms:
        molecule.AddAtom(Chem.Atom(atom.element.symbol))

    for atom_1, atom_2 in topology.bonds:
        molecule.AddBond(atom_1.index, atom_2.index)

    return molecule

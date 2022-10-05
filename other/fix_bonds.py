from rdkit.Chem import AllChem


def fix_bond_order_from_smiles(molecule, smiles):
    """ Assign the correct bond order to a molecule with the corresponding smiles.

        Parameters
        ----------
        molecule : rdkit.Mol
            A molecule with incorrect bond orders.

        smiles: str
            The smiles that is used as a template to fix bond orders

        Returns
        -------
        rdkit.Mol
            Molecule with correct bond orders
    """
    template = AllChem.MolFromSmiles(smiles)
    return AllChem.AssignBondOrdersFromTemplate(template, molecule)

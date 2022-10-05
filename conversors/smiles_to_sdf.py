from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_sdf(file_name, smiles):
    """ Creates a sdf file from a list of smiles. """
    molecules = []
    for code in smiles:
        mol = Chem.MolFromSmiles(code)
        assert mol is not None, "Failed to convert smiles to mol"
        mol = Chem.AddHs(mol)
        AllChem.EmbedMultipleConfs(mol, 1)
        molecules.append(mol)

    writer = Chem.SDWriter(file_name)
    for mol in molecules:
        writer.write(mol)

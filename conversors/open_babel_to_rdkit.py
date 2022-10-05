from io import StringIO, BytesIO
from rdkit import Chem


def open_babel_to_rdkit(molecule):
    """ Transform an openbabel molecule to rdkit. """
    # TODO: transform without writing a file
    sdf_str = molecule.write("sdf")
    ligand_sio = StringIO(sdf_str)
    ligand_bio = BytesIO(ligand_sio.read().encode("utf8"))
    ligand = [mol for mol in Chem.ForwardSDMolSupplier(ligand_bio)][0]
    return ligand

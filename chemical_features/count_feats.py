from chemical_features import find_smarts_matches, get_chemical_features, smarts_patterns
from rdkit import Chem
from pprint import pprint


def count_chemical_features(file_name):

    molecule = Chem.MolFromPDBFile(file_name)
    matches = find_smarts_matches(molecule, smarts_patterns)
    count = {k: len(v) for k, v in matches.items()}
    print("Chemical features found with smarts pattern matching:")
    pprint(count)

    feats = get_chemical_features(molecule)
    count = {k: len(v) for k, v in feats.items()}
    print("Chemical features found with rdkit feature factory:")
    pprint(count)


if __name__ == "__main__":
    test_file_1 = "../test_cases/eralpha/1qku/1qku.pdb"
    test_file_2 = "../test_cases/rhinovirus/1ncr/1ncr.pdb"
    count_chemical_features(test_file_1)
    count_chemical_features(test_file_2)
    print("DONE")

from rdkit import RDConfig, Chem
from rdkit.Chem import ChemicalFeatures
import os
import numpy as np
from collections import defaultdict


smarts_patterns = {
    '*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]': 'Hydrophobe',
    'N=[CX3](N)-N': 'PosIonizable',
    '[#16!H0]': 'Donor',
    '[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]': 'Donor',
    '[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]': 'Acceptor',
    '[#8!H0&!$([OH][C,S,P]=O)]': 'Donor',
    '[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,'
    'CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,'
    'I]': 'Hydrophobe',
    '[$([+,+2,+3])&!$(*[-,-2,-3])]': 'PosIonizable',
    '[$([-,-2,-3])&!$(*[+,+2,+3])]': 'NegIonizable',
    '[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]': 'Hydrophobe',
    '[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,'
    'CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]': 'Hydrophobe',
    '[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]': 'Hydrophobe',
    '[$([CX3,SX3,PX3](=O)[O-,OH])](=O)[O-,OH]': 'NegIonizable',
    '[$([CX3](=N)(-N)[!N])](=N)-N': 'PosIonizable',
    '[$([NX3]([CX4])([CX4,#1])[CX4,#1])&!$([NX3]-*=[!#6])]': 'PosIonizable',
    '[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]': 'Acceptor',
    '[$([SX4,PX4](=O)(=O)[O-,OH])](=O)(=O)[O-,OH]': 'NegIonizable',
    '[$([S]~[#6])&!$(S~[!#6])]': 'Hydrophobe',
    '[C&r3]1~[C&r3]~[C&r3]1': 'Hydrophobe',
    '[C&r4]1~[C&r4]~[C&r4]~[C&r4]1': 'Hydrophobe',
    '[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1': 'Hydrophobe',
    '[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1': 'Hydrophobe',
    '[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1': 'Hydrophobe',
    '[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1': 'Hydrophobe',
    '[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]': 'Hydrophobe',
    'a1aaaa1': 'Aromatic',
    'a1aaaaa1': 'Aromatic',
    'c1nn[nH1]n1': 'NegIonizable'}


def load_smarts_fdef(file_name: str):
    """ Load custom chemical feature definitions from a txt file.
    
        Feature definitions are SMART strings with their respective chemical feature
        name. 

        Parameters
        ----------
        file_name : str
            Name of the file containing the smarts feature definitions
        
        Returns
        -------
        features : dict
            Dictionary which keys are SMARTS strings and values are feature names

        Note
        -----
        The file must contain a SMARTS string follwed by the feature name. 
        
        Example:

        # Aromatic
        a1aaaaa1 Aromatic
        a1aaaa1 Aromatic

        Lines started with # are considered as comments
    """
    features = {}
    # Load custom features file
    with open(file_name, "r") as file:
        for line in file:
            if line[0] == "#" or line[0] == '\n':
                continue
            line = line.split(' ')
            feat_def = line[0]  # The smarts string
            feat_name = line[1].rstrip()
            features[feat_def] = feat_name

    return features


def find_smarts_matches(molecule, smarts_patterns):
    """ Find matches of smarts patterns in a molecule

        Parameters
        ----------

        Returns
        -------


    """
    features = defaultdict(list)

    for smarts, feat_name in smarts_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        atom_indices = molecule.GetSubstructMatches(pattern)
        if len(atom_indices) > 0:
            for indices in atom_indices:
                features[feat_name].append(indices)

    return features


def get_chemical_features(molecule):
    """ Find chemical features in a molecule using rdkit feature factory.
    """
    features = defaultdict(list)
    # Use rdkit feature factory
    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    featdef = ChemicalFeatures.BuildFeatureFactory(fdefName)
    chemical_features = featdef.GetFeaturesForMol(molecule)
    for feature in chemical_features:
        feat_name = feature.GetFamily()
        atom_indices = feature.GetAtomIds()
        if len(atom_indices) > 0:
            for indices in atom_indices:
                features[feat_name].append(indices)

    return features


def get_feature_coords_and_direction(ligand, atom_indices,
                                     conformer_index, directionality):
    """ Obtain the coordinates and if specified the directionality of a chemical feature.

        Parameters
        ----------
        ligand : rdkit.Chem.Mol
            A ligand

        conformer_index : int
            The conformer whose coordinates will be used to obtain the pharmacophoric
            points.

        directionality : bool
            Whether to compute the direction vector of that point.

        Returns
        -------


    """
    if len(atom_indices) > 1:
        # Find the centroid
        # Aromatic, hydrophobic, positive or negative feature
        coords = _feature_centroid(ligand, atom_indices, conformer_index)
        # Find direction vector
        if directionality:
            direction = _aromatic_direction_vector(ligand, atom_indices,
                                                   conformer_index)
        else:
            direction = None
    else:
        # Find the centroid
        # Donor or acceptor feature
        position = ligand.GetConformer(conformer_index).GetAtomPosition(atom_indices[0])
        coords = np.zeros((3,))
        coords[0] = position.x
        coords[1] = position.y
        coords[2] = position.z
        # Find direction vector
        if directionality:
            direction = _donor_acceptor_direction_vector(ligand, atom_indices[0],
                                                         coords, conformer_index)
        else:
            direction = None

    return coords, direction


def _feature_centroid(molecule, atom_indxs, conformer_index):
    """
        Get the 3D coordinates of the centroid of a feature that encompasses more than
        one atom. This could be aromatic, hydrophobic, negative and positive features

        Parameters
        ----------
        molecule : rdkit.Chem.Mol
                Molecule that contains the feature which centroid will be computed

        atom_indxs : tuple of int
                Indices of the atoms that belong to the feature

        conformer_index : int
                Index of the conformer for which the feature centroid will be computed

        Returns
        -------
        centroid : numpy.ndarray of shape (3, )
            Array with the coordinates of the centroid of the feature.

    """

    n_atoms = len(atom_indxs)
    coords = np.zeros((n_atoms, 3))
    for j, idx in enumerate(atom_indxs):
        position = molecule.GetConformer(conformer_index).GetAtomPosition(idx)
        coords[j, 0] = position.x
        coords[j, 1] = position.y
        coords[j, 2] = position.z

    centroid = coords.mean(axis=0)

    return centroid


def _donor_acceptor_direction_vector(molecule: Chem.Mol, feat_type: str, atom_indx: int,
                                     coords: np.ndarray, conformer_idx: int) -> np.ndarray:
    """
        Compute the direction vector for an H bond donor or H bond acceptor feature

        Parameters
        ----------
        molecule : rdkit.Mol
                Molecule that contains the feature which direction vector will be computed.

        feat_type : str
                Type of feature. Whether is a donor or acceptor.

        atom_indx : int
                Index of the H bond acceptor or donor atom.

        coords : numpy.ndarray; shape(3,)
                Coordinates of the H bond acceptor or donor atom.

        conformer_idx : int
                Index of the conformer for which the direction vector will be computed.

        Returns
        -------
        direction : numpy.ndarray; shape(3,)
                Coordinates of the direction vector.

    """
    direction = np.zeros((3,))
    atom = molecule.GetAtomWithIdx(atom_indx)
    for a in atom.GetNeighbors():
        if a.GetSymbol() == "H":
            continue
        position = molecule.GetConformer(conformer_idx).GetAtomPosition(a.GetIdx())
        direction[0] += position.x - coords[0]
        direction[1] += position.y - coords[1]
        direction[2] += position.z - coords[2]
    if feat_type == "Donor":
        direction = -direction
    return direction


def _aromatic_direction_vector(molecule, atom_indxs, conformer_idx):
    """ Compute the direction vector for an aromatic feature.

        Parameters
        ----------
        molecule : rdkit.Mol
                Molecule that contains the feature which direction vector will be computed.

        atom_indxs : tuple of int
                Indices of the aromatic atoms.

        conformer_idx : int
                Index of the conformer for which the direction vector will be computed.

        Returns
        -------
        direction : numpy.ndarray; shape(3,)
                Coordinates of the direction vector.

    """
    coords = np.zeros((3, 3))  # Take just the first three atoms
    for j, idx in enumerate(atom_indxs[0:3]):
        position = molecule.GetConformer(conformer_idx).GetAtomPosition(idx)
        coords[j, 0] = position.x
        coords[j, 1] = position.y
        coords[j, 2] = position.z

    # Find the vector normal to the plane defined by the three atoms
    u = coords[1, :] - coords[0, :]
    v = coords[2, :] - coords[0, :]
    direction = np.cross(u, v)

    return direction


def add_features_to_view(view, molecule, features):
    """ Adds the chemical features to a ngl view. They are represented as spheres."""
    palette = {
        'PosIonizable': '#3498DB',  # Blue
        'NegIonizable': '#884EA0',  # Purple
        'Acceptor': '#B03A2E',  # Red
        'Donor': '#17A589',  # Green
        'Hydrophobe': '#F5B041',  # Orange
        'Aromatic': '#F1C40F',  # Yellow
    }

    for feat_name, atom_indices in features.items():

        try:
            color = palette[feat_name]
        except KeyError:
            continue

        for indices in atom_indices:
            centroid = _feature_centroid(molecule, indices, 0)
            n_components = len(view._ngl_component_ids)
            view.shape.add_sphere(centroid, color, 0.5, feat_name)
            view.update_representation(component=n_components, repr_index=0, opacity=0.5)

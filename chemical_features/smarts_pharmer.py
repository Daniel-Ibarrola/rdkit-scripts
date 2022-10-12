from pprint import pprint
# Smarts patterns definitions. Taken from Pharmer

# The following are for ligands

aromatic = [
    "a1aaaaa1",
    "a1aaaa1",
]

hydrogen_donor = [
    "[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]",
    "[#8!H0&!$([OH][C,S,P]=O)]",
    "[#16!H0]",
]

hydrogen_acceptor = [
    "[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]",
    "[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]",
]

positive_ion = [
    "[+,+2,+3,+4]",
    # amidine
    "[$(CC)](=N)N",
    # guanidine
    "[$(C(N)(N)=N)]",
    "[$(n1cc[nH]c1)]",
]

negative_ion = [
    "[-,-2,-3,-4]",
    "C(=O)[O-,OH,OX1]",
    "[$([S,P](=O)[O-,OH,OX1])]",
    "c1[nH1]nnn1",
    "c1nn[nH1]n1",
    "C(=O)N[OH1,O-,OX1]",
    "C(=O)N[OH1,O-]",
    "CO(=N[OH1,O-])",
    "[$(N-[SX4](=O)(=O)[CX4](F)(F)F)]",
]

hydrophobic = [
    # branched terminals as one point
    "[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
    "[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,"
    "CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,"
    "I]",
    "*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
    # simple rings only; need to combine points to get good results for 3d structures
    "[C&r3]1~[C&r3]~[C&r3]1",
    "[C&r4]1~[C&r4]~[C&r4]~[C&r4]1",
    "[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1",
    "[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1",
    "[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1",
    "[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1",
    # aliphatic chains
    "[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
    "[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
    "[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,"
    "CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
    # sulfur (apparently)
    "[$([S]~[#6])&!$(S~[!#6])]",
]

# The following are for proteins

positive_ion_protein = [
    "[+,+2,+3,+4]",
    # amidine
    # guanidine
    "[$(C(N)(N)=N)]", "[$(n1cc[nH]c1)]",
]

negative_ion_protein = [
    "[-,-2,-3,-4]",
    "C(=O)[O-,OH,OX1]",
]

hydrophobic_protein = [
    # branched terminals as one point
    "[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
    "[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,"
    "CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,"
    "I]",
    "[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
    "[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
    "[$([S]~[#6])&!$(S~[!#6])]",
]


def create_dict_ligand_smarts():

    ligand_smarts = {}

    for smarts in aromatic:
        ligand_smarts[smarts] = "Aromatic"
    for smarts in hydrogen_donor:
        ligand_smarts[smarts] = "Donor"
    for smarts in hydrogen_acceptor:
        ligand_smarts[smarts] = "Acceptor"
    for smarts in positive_ion:
        ligand_smarts[smarts] = "PosIonizable"
    for smarts in negative_ion:
        ligand_smarts[smarts] = "NegIonizable"
    for smarts in hydrophobic:
        ligand_smarts[smarts] = "Hydrophobe"

    return ligand_smarts


def create_dict_protein_smarts():

    protein_smarts = {}
    for smarts in aromatic:
        protein_smarts[smarts] = "Aromatic"
    for smarts in hydrogen_donor:
        protein_smarts[smarts] = "Donor"
    for smarts in hydrogen_acceptor:
        protein_smarts[smarts] = "Acceptor"
    for smarts in positive_ion_protein:
        protein_smarts[smarts] = "PosIonizable"
    for smarts in negative_ion_protein:
        protein_smarts[smarts] = "NegIonizable"
    for smarts in hydrophobic_protein:
        protein_smarts[smarts] = "Hydrophobe"

    return protein_smarts


pprint(create_dict_protein_smarts())

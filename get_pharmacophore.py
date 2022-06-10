from pmapper.pharmacophore import Pharmacophore as P
from rdkit import Chem
from rdkit import RDConfig
import os
import pandas as pd
from Bio.SeqUtils import seq3

# Иницилизируем фармакофоры, которые будут использоваться в RDkit, возможно файл BaseFeatures.fdef можно отредактировать, чтобы там были нужные нам фармакофоры, но у меня не получилось

fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef') 
factory = Chem.ChemicalFeatures.BuildFeatureFactory(fdefName)

# SMILES-представление для аминокислот
aa_smiles = {'ALA': 'O=C(O)[C@H](C)N',
             'CYS': 'O=C(O)[C@H](CS)N', 
             'ASP': 'O=C(O)[C@H](CC(O)=O)N', 
             'GLU': 'O=C(O)[C@H](CCC(O)=O)N', 
             'PHE': 'O=C(O)[C@H](CC1=CC=CC=C1)N', 
             'GLY': 'O=C(O)CN', 
             'HIS': 'O=C(O)[C@H](CC1=CNC=N1)N', 
             'ILE': 'O=C(O)[C@H]([C@@H](C)CC)N', 
             'LYS': 'O=C(O)[C@H](CCCCN)N', 
             'LEU': 'O=C(O)[C@H](CC(C)C)N', 
             'MET': 'O=C(O)[C@H](CCSC)N', 
             'ASN': 'O=C(O)[C@H](CC(N)=O)N', 
             'PRO': 'O=C(O)[C@H]1NCCC1', 
             'GLN': 'O=C(O)[C@H](CCC(N)=O)N', 
             'ARG': 'O=C(O)[C@H](CCCNC(N)=N)N', 
             'SER': 'O=C(O)[C@H](CO)N', 
             'THR': 'O=C(O)[C@H]([C@H](O)C)N', 
             'VAL': 'O=C(O)[C@H](C(C)C)N', 
             'TRP': 'O=C(O)[C@H](CC1=CNC2=CC=CC=C12)N',
             'TYR': 'O=C(O)[C@H](CC1=CC=C(C=C1)O)N'}

# Сопоставление pmapper фармакофоров и rdkit фармакофоров
compare = {"A": "Acceptor",
           "D": "Donor",
           "P": "PosIonizable",
           "N": "NegIonizable",
           "H": "Hydrophobe",
           "a": "Aromatic"}

# Считаем фармакофоры для всех аминокислот с помощью двух тулов. Сохраняем результаты в 2 датафрейма

pharmacophore_RDkit = pd.DataFrame(columns = [aa for aa in aa_smiles], index = ["Hydrophobe", "PosIonizable", "NegIonizable", "Acceptor", "Donor", "Aromatic", "Sulphur"])
pharmacophore_RDkit = pharmacophore_RDkit.fillna(0)
pharmacophore_Pmapper = pd.DataFrame(columns = [aa for aa in aa_smiles], index = ["Hydrophobe", "PosIonizable", "NegIonizable", "Acceptor", "Donor", "Aromatic", "Sulphur"])
pharmacophore_Pmapper = pharmacophore_Pmapper.fillna(0)
for aa in aa_smiles:
    #RDkit
    mol = Chem.MolFromSmiles(aa_smiles[aa])
    mol = Chem.AddHs(mol) 
    Chem.AllChem.EmbedMolecule(mol, randomSeed=42)
    feats = factory.GetFeaturesForMol(mol)
    for i in range(len(feats)):            
        ph = feats[i].GetFamily()
        if feats[i].GetFamily() == "LumpedHydrophobe":
            ph = "Hydrophobe"
        pharmacophore_RDkit[aa][ph] += 1
    #Pmapper
    p = P()
    p.load_from_mol(mol)
    pharmacophore_dict = p.get_features_count()
    for i in pharmacophore_dict:
        pharmacophore_Pmapper[aa][compare[i]] += pharmacophore_dict[i]
    if "S" in aa_smiles[aa]:
        pharmacophore_RDkit[aa]["Sulphur"] += 1
        pharmacophore_Pmapper[aa]["Sulphur"] += 1
        
# Функции для подсчета разницы фармакофоров        
def get_pharm_diff_RDkit(row):
    wild = seq3(row["WILD_TYPE"]).upper()
    mutant = seq3(row["MUTANT"]).upper()
    return list(pharmacophore_RDkit[mutant] - pharmacophore_RDkit[wild])
def get_pharm_diff_Pmapper(row):
    wild = seq3(row['WILD_TYPE']).upper()
    mutant = seq3(row['MUTANT']).upper()
    return list(pharmacophore_Pmapper[mutant] - pharmacophore_Pmapper[wild])

# Загружаем наш датасет с посчитанными CSM после create_CSM.py
PDB_dataset = pd.read_csv("/path/to/dataset_with_CSM.csv")

# Считаем фармакофоры двумя способами

PDB_dataset["pharmacophore_RDkit"] = PDB_dataset.apply(get_pharm_diff_RDkit, axis=1)
PDB_dataset["pharmacophore_Pmapper"] = PDB_dataset.apply(get_pharm_diff_Pmapper, axis=1)

#Сохраняем датасет с CSM и факрмакофорами

PDB_dataset.to_csv("/path/to/folder/dataset_with_CSM_and_PH.csv")

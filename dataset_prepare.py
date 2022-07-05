import csv
import pandas as pd
import numpy as np
import glob
from Bio.PDB import PDBParser, PDBIO, Select # for needed model and chain
from pdbtools import pdb_fetch # for raw PDB download
import subprocess # for run bash script
from schrodinger import structure # for PRR
from schrodinger.protein import reliability # for PRR

# Upload dataset from CSV to pd.DataFrame 

PDB_dataset = pd.read_csv("/path/to/datadrame.csv")

'''We should have such columns 
'PDB' (PDB_id; str), 
'WILD_TYPE' (naive one-letter residue; str), 
'CHAIN' (chain, where mutation occurs; str),
'POSITION' (mutation residue number, could be with letters; str), 
'MUTANT' (mutant one-letter residue; str), 
'DDG' (experimental ddG; float)
If columns are different - delete unneeded and rename left as above'''

# Drop duplicates and save only mean value of ddG for one mutation (if there are values with different ddG sign for one mutation - delete it)

PDB_dataset = PDB_dataset.drop_duplicates()
PDB_dataset = PDB_dataset.groupby(['PDB', 'CHAIN', 'WILD_TYPE', 'POSITION', 'MUTANT'], as_index=False).filter(lambda x:(x.ddG.all() > 0) || (x.ddG.all() < 0))
PDB_dataset = PDB_dataset.groupby(['PDB', 'CHAIN', 'WILD_TYPE', 'POSITION', 'MUTANT'], as_index=False).aggregate({"ddG":"mean"})



# Create extra column PDB_CHAIN, will need it later

PDB_dataset["PDB_CHAIN"] = PDB_dataset["PDB"] + PDB_dataset["CHAIN']

# Get unique PDB_id from dataset 

PDB_unique = PDB_dataset["PDB"].unique()

# Download raw PDB 

for PDB_id in PDB_unique:
    PDB_file = pdb_fetch.run(PDB_id[:4])
    with open(f"/path/to/folder/{PDB_id[:4]}.pdb",mode="a") as f:
        for i in PDB_file:
            f.write(i)
            
# Get unique PDB_CHAIN from dataset 

PDB_id_chain_unique = PDB_dataset['PDB_CHAIN'].unique()

                                                               
# Create class to remove heteroatoms and non-protein atoms   
                                                               
class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0                         
                                                               
# Save only 1 model and needed chain (or chains if several for one protein)

parser = PDBParser()
io = PDBIO()
for PDB_id_chain in PDB_id_chain_unique:    
    structure = parser.get_structure(PDB_id_chain[:4], f"/path/to/folder/{PDB_id_chain[:4]}.pdb")
    chain = structure[0][PDB_id_chain[4]] # get 1 model and needed chain
    io.set_structure(chain)
    io.save(f"/path/to/folder/{PDB_id_chain[:5]}_model_1.pdb", NonHetSelect())

# Apply Protein Preparation Wizard (PPW) to all structures (can be done only in bash). Make sure you can run Schrodinger from command line
# As a result we get {PDB_CHAIN}_prep.pdb files
                                                               
import subprocess
subprocess.run('dataset_prepwiz.sh')

# Check structures after PPW by Protein Reliability Report (PRR) for missing atoms, missung chains and steric clashes

errors = {}
for pdb_file in glob.glob("*_prep.pdb"):
    st = structure.StructureReader.read(pdb_file)
    check = reliability.run_quick_check(st)
    if check.getTextSummary() != []:
        errors[pdb_file] = check.getTextSummary()
structures_with_errors = []
for structure, report in errors.items():
    problems = ['missing loops','missing atoms','steric clashes']
    matching = [report_item for report_item in report if any(problem in report_item for problem in problems)]
    if matching != []:
        print("There is a problem with this structure:", structure)
        structures_with_errors.append(structure[:5])

# Delete all structures from previous step from dataset
                                                            
for structure in structures_with_errors:
    PDB_dataset.drop(PDB_dataset[PDB_dataset['PDB_CHAIN'] == structure].index, inplace=True)

# Save filtered dataset
                                                            
PDB_dataset.to_csv(f'/path/to/folder/dataset.csv', header=True, index=False, mode='w')

# As a result we get dataset.csv file with needed columns and bunch of {PDB_CHAIN}_prep.pdb files

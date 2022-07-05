import glob
import dhutil
import pandas as pd
import csv
from schrodinger import structure
from Bio.SeqUtils import seq3

# .pdb to .mae

for pdb_file in glob.glob("*_prep.pdb"):
    st = structure.StructureReader.read(pdb_file)
    modified_fname = pdb_file[0:5] + '.mae'
    with structure.StructureWriter(modified_fname) as writer:
        writer.append(st)

# Create folder for each PDB_CHAIN, create .txt file with list of mutations and .sh file with command. Transfer .mae, .txt and .sh to relevant folder

PDB_dataset = pd.read_csv('/path/to/dataset.csv')
for PDB_id_chain in PDB_dataset.PDB_CHAIN.unique():
    os.mkdir(PDB_id_chain[0:5])
    shutil.copy(f'/path/to/{PDB_id_chain[0:5]}.mae', f'/path/to/{PDB_id_chain[0:5]}/{PDB_id_chain[0:5]}.mae')
    temp_dataset = PDB_dataset.loc[(PDB_dataset['PDB_CHAIN'] == PDB_id_chain)]
    temp_dataset["to_save"] = PDB_dataset['CHAIN'] + ":" + PDB_dataset['POSITION'] + " " + PDB_dataset['MUTANT'].map(seq3).str.upper()
    temp_dataset["to_save"].to_csv(f"/path/to/{PDB_id_chain[0:5]}/{PDB_id_chain[0:5]}.txt" , header=False, index=None, quoting=csv.QUOTE_NONE, mode='w', quotechar='', escapechar = "\\")
    mutations_filename = PDB_id_chain[0:5] + ".txt"
    residue_scanning_file = PDB_id_chain[0:5] + ".mae"
    command = '"${SCHRODINGER}/run" residue_scanning_backend.py -WAIT -jobname residue_scanning_CSM -fast -res_file ' + mutations_filename + " -refine_mut prime_residue -calc hydropathy,rotatable,vdw_surf_comp,sasa_polar,sasa_nonpolar,sasa_total -dist 0.00 " + residue_scanning_file + " -add_res_scan_wam -HOST localhost:1 -NJOBS 1 -TMPLAUNCHDIR"
    with open(f'/path/to/{PDB_id_chain[0:5]}/{PDB_id_chain[0:5]}.sh', "w") as text_file:
        print(command, file=text_file)

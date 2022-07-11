import glob
import dhutil
import pandas as pd
import csv
import re
import os
from schrodinger import structure
from Bio.SeqUtils import seq3, seq1

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

# After running RS, put you result folders named PDB_CHAIN into RS_done folder. Now we add DDG_RS column to dataset with ddG values after RS

os.chdir(f'/path/to/RS_done')
PDB_dataset["DDG_RS"] = 0
for pdb_file in glob.glob("*"):
    print(pdb_file)    
    os.chdir(f'/path/to/RS_done/{pdb_file}')
    PDB_dataset_temp = PDB_dataset.loc[PDB_dataset["PDB_CHAIN"] == pdb_file]
    for ind in number_of_mut:
        st = structure.StructureReader.read('residue_scanning_CSM-out.maegz', index=ind+2)
        mutation = st.property['s_bioluminate_Mutations']
        mutation_split = re.search('(.):(.*)\((.*)->(.*)\)', mutation)
        chain = mutation_split.group(1)
        position = mutation_split.group(2)
        wt = seq1(mutation_split.group(3), undef_code='H')
        mut = seq1(mutation_split.group(4), undef_code='H')
        index = PDB_dataset.loc[(PDB_dataset["CHAIN"] == chain) &
                        (PDB_dataset["POSITION"] == position) &
                        (PDB_dataset["WILD_TYPE"] == wt) &
                        (PDB_dataset["MUTANT"] == mut)].index
        PDB_dataset.at[index,"DDG_RS"] = -st.property['r_psp_Prime_delta_Stability']
    os.chdir(f'/Users/georgii_mamistvalov/Documents/Projects/CSM/Datasets_stability/S2648_after_prepwizard_done')

# Save updated dataset

PDB_dataset.to_csv(f'/path/to/dataset.csv', header=True, index=False, mode='w')

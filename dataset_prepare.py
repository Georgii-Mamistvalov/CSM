# Цель данного блока - получить набор pdb-структур из датасета (в фомате .csv). Также мы переводим изначальный датасет в общий формат (единые столбцы).

import csv
import os
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser # для выделения нужной цепи и нужной модели
from Bio.PDB import PDBIO # для сохранения результата
from Bio.PDB import Select # для выбора только нужной цепи
from pdbtools import pdb_fetch # для скачивания PDB

# Загружаем датасет из CSV в pd.DataFrame 

PDB_dataset = pd.read_csv("/path/to/datadrame.csv")

'''Предполагается наличие столбцов типа 
'PDB' (PDB_id; str), 
'WILD_TYPE' (нативный остаток однобуквенный; str), 
'CHAIN' (цепь, в которой выполняется замена; str),
'POSITION' (номер заменяемого остатка; int), 
'MUTANT' (мутантный остаток однобуквенный; str), 
'DDG' (экспериментально полученное ddG, с ним будем сравнивать и на нем будем учить модель; float)
Если столбцы не такие - стоит удалить ненужные и переименовать нужные в указанные выше'''

# Создаем доп. столбец PDB+CHAIN, понадобится для этапа (# Обрабатываем, сохраняя только первую модель и нужную цепь)

PDB_dataset["PDB_id_chain"] = PDB_dataset["PDB"] + PDB_dataset["CHAIN']

# Получаем список уникальных PDB_id из датасета 

PDB_unique = PDB_dataset.PDB.unique()

# Скачиваем все уникальные PDB в первоначальном виде

for PDB_id in PDB_unique:
    PDB_file = pdb_fetch.run(PDB_id[0:4])
    with open("/path/to/folder/"+PDB_id[0:4] + '.pdb',mode="a") as f:
        for i in PDB_file:
            f.write(i)
            
# Получаем список уникальных PDB_id_chain из датасета 

PDB_id_chain_unique = PDB_dataset.PDB_id_chain.unique()

                                                               
# Делаем класс для удаления гетероатомов и всего, кроме белка (ВРОДЕ ТАК, НАДО ЭТО УТОЧНИТЬ)  
                                                               
class NonHetSelect(Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0                         
                                                               
# Обрабатываем, сохраняя только первую модель и нужную цепь.

parser = PDBParser()
io = PDBIO()
for PDB_id_chain in PDB_id_chain_unique:    
    structure = parser.get_structure(PDB_id_chain[0:4], "/path/to/folder/"+PDB_id_chain[0:4] + '.pdb')
    chain = structure[0][PDB_id_chain[4]] # берем первую модель и нужную цепь
    io.set_structure(chain)
    io.save("/path/to/folder/"+PDB_id_chain[0:5]+"_model_1.pdb", NonHetSelect())

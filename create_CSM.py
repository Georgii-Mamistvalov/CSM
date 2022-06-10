import csv
import numpy as np
import pandas as pd
from graphein.protein.graphs import read_pdb_to_dataframe, process_dataframe, initialise_graph_with_metadata, add_nodes_to_graph
from graphein.protein.subgraphs import extract_subgraph_from_point
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
from graphein.protein.edges.atomic import add_atomic_edges
import re
from Bio.SeqUtils import seq3
import math

# функция для подсчета матрицы попарных расстояний между вершинами графа g

def compute_length_of_edges(g):
    length = []
    for a1 in g.nodes:
        for a2 in g.nodes:
            if a1 != a2:
                dist = math.dist(g.nodes[a1]["coords"], g.nodes[a2]["coords"])
                length.append(dist)
    return length

# функция для подсчета вектора CSM (распределение расстояний), lenght - матрица попарных расстояний между атомами, dmin - начальное расстояние,
# dmax - максимальное расстояние, dstep - шаг

def compute_CSM(lenth, dmin, dmax, dstep):
    CSM = list(0 for i in np.arange(dmin, dmax, dstep))
    i = 0
    d = dmin
    while d <= dmax - dstep:
        for dist in lenth:
            if dist >= d and dist <= d+dstep:
                CSM[i] += 1
        d += dstep
        i += 1
    return CSM

# основная функция - представляет белок в виде графа и использует compute_length_of_edges и compute_CSM для подсчета CSM
def get_final_CSM(row):
    PDB_id = row["PDB"]
    chain = row["CHAIN"]
    aminoacid = row["WILD_TYPE"]
    position = row["POSITION"]
    
    params_to_change = {"granularity": "atom", "edge_construction_functions": [add_atomic_edges], 'insertions': True} # Устанавливаем необходимые параметры для атомарного графа
    config = ProteinGraphConfig(**params_to_change) # меняем параметры в конфиге построения графа
    
    pdb_path ="/path/to/folder/"+PDB_id+chain+"_model_1.pdb" # путь до pdb
    
    raw_df = read_pdb_to_dataframe(pdb_path=pdb_path, granularity = "atom") # представляем pdb в виде pd.DataFrame. Тут можно при желании выбрать, в каком виде хотим представить белок - атомы, Сa и тд.
    df = process_dataframe(raw_df, insertions = True, granularity = "atom") # препроцессим (ВАЖНО на этом этапе также можно задать granularity, надо это потестить, пока просто делаем аналогично строке выше)
    
    df["node_id"] = [number + "_" + letter for number, letter in zip(df["node_id"], df["insertion"])] # меняем названия вершин, иначе 27А и 27С имеют одинаковое node_id и часть вершин может не нарисоваться
    
    g = initialise_graph_with_metadata(protein_df=df, raw_pdb_df=raw_df.df["ATOM"], pdb_id = PDB_id, granularity = "atom") #иницилизируем граф 
    g = add_nodes_to_graph(g) # добавляем вершины
    
    mutation_ID = chain+":"+seq3(aminoacid).upper()+":"+ re.sub(r'[^\d]', '', position) + ":CA" +  "_" + re.sub(r'[^a-zA-Z]', '', position) # делаем node_id нашей замены в нужном виде, считаем атомы вокруг Ca мутируемой аминокислоты
    #считаем координаты
    mutation_point_x = g.nodes[mutation_ID]["coords"][0]
    mutation_point_y = g.nodes[mutation_ID]["coords"][1]
    mutation_point_z = g.nodes[mutation_ID]["coords"][2]
    
    # задаем параметры подсчета CSM (в ангстремах)
    d_max = 10
    d_step = 0.2
    
    # выделяем подграф вокруг нашей мутации радиуса d_max (wild-type residue environment)
    s_g = extract_subgraph_from_point(g, centre_point=(mutation_point_x, mutation_point_y, mutation_point_z), radius=d_max)
    
    # считаем CSM для нашего подграфа
    length = compute_length_of_edges(s_g) # вектор с попарными расстояниями между атомами. Матрица тут не нужна, тк считаем потом просто количества разных расстояний
    CSM = compute_CSM(length, 0, d_max, d_step) # вектор количества расстояний в выбранных диапазонах
    CSM_pdf = CSM / sum(CSM) # распределение вероятности (Probability Distribution Function)
    CSM_cdf = np.cumsum(CSM_pdf) # кумулятивное распределение вероятности (Сumulative Distribution Function)
    
    # Есть вопрос, что же возвращать - CSM, CSM_pdf или CSM_cdf. Авторы используют CSM_cdf. Мне кажется более интересно использовать CSM или СSM_pdf 
    # НАДО ТЕСТИТЬ
    
    return CSM_cdf

# Подгружаем датасет

PDB_dataset = pd.read_csv("/path/to/dataset.csv")

# Добавляем CSM в наш датасет

PDB_dataset["CSM"] = PDB_dataset.apply(get_final_CSM, axis=1)

# Сохраняем обновленный датасет

PDB_dataset.to_csv(r'/path/to/folder/dataset_with_CSM_name.csv', header=True, mode='w')

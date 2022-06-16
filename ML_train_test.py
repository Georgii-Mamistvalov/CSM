from sklearn import tree
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.metrics import mean_squared_error
from sklearn.feature_selection import r_regression
from sklearn.ensemble import ExtraTreesRegressor

sns.set_theme(color_codes=True)

PDB_dataset = pd.read_csv("/path/to/folder/dataset_with_CSM_and_PH.csv")
'''
'PDB' (PDB_id; str), 
'WILD_TYPE' (нативный остаток однобуквенный; str), 
'CHAIN' (цепь, в которой выполняется замена; str),
'POSITION' (номер заменяемого остатка; int), 
'MUTANT' (мутантный остаток однобуквенный; str), 
'DDG' (экспериментально полученное ddG, с ним будем сравнивать и на нем будем учить модель; float)
'PDB_id_chain' (PDB+CHAIN, str)
'DDG_web_tool' (предсказанное ddG веб-тулом авторов статьи; float)
'CSM' (вектор-список CSM, str)
'pharmacophore_RDkit' (вектор-список фармакофоров посчитанных по RDkit, str)
'pharmacophore_Pmapper' (вектор-список фармакофоров посчитанных по Pmapper, str)
'''

# Конвертируем str в list
df["pharmacophore_RDkit"] = df["pharmacophore_RDkit"].map(json.loads)
df["pharmacophore_Pmapper"] = df["pharmacophore_Pmapper"].map(json.loads)
df["CSM"] = df["CSM"].map(json.loads)

# Готовим список названий колонок для дальнейшего разделения списков на столбцы (каждое значение свой столбец)
list_names_CSM = []
for i in range(len(df["CSM"][0])):
    list_names_CSM.append("CSM_"+str(i))
list_names_pharmacophore_RDkit = []
for i in range(len(df.pharmacophore_RDkit[0])):
    list_names_pharmacophore_RDkit.append("pharmacophore_RDkit_"+str(i))
    
# Готовим датафрейм X с количеством столбцов n + m (n - длина вектора CSM, m - длина вектора Pharmacophore) и вектор y с результатами эксперимента и результатами тула
X = df[["CSM","pharmacophore_RDkit"]]
y = df[["DDG", "DDG_web_tool"]]
X[list_names_CSM] = pd.DataFrame(X.CSM.tolist(), index= X.index)
X[list_names_pharmacophore_RDkit] = pd.DataFrame(X.pharmacophore_RDkit.tolist(), index= X.index)
X.head()
X = X.drop(['CSM',"pharmacophore_RDkit"], axis=1)

# Создаем и тренируем дерево. Ищем оптимальные параметры с помощью RandomizedSearchCV

etr = ExtraTreesRegressor(random_state = 0)
X_train, X_test, y_train, y_test = train_test_split(X, y["DDG"], random_state=0)
params = {
    'max_depth': [i for i in range(1, 51, 1)],
    'min_samples_split': [i for i in range(2, 10, 2)],
    'min_samples_leaf':[i for i in range(1, 8)],
    'n_estimators':[i for i in range(10, 301, 10)]
}
search = RandomizedSearchCV(etr, cv=20, param_distributions=params, n_jobs=-1, n_iter=100)
search.fit(X_train, y_train)

etr = search.best_estimator_

prediction = etr.predict(X_test)

ddG_all = {"DDG_pred_home": prediction, "DDG_exp": y_test["DDG"], "DDG_pred_web": y_test["DDG_web_tool"]}



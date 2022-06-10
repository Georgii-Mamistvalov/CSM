
# Цель - получить значения ddG для исследуемого датасета из веб-инструмента, реализованного авторами статьи. 

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from bs4 import BeautifulSoup
from webdriver_manager.chrome import ChromeDriverManager
import time
from io import StringIO

# Для начала необходимо понять, какие фаилы веб-инструмент принимает на вход. В случае со стабильностью - pdb-файл + список мутаций.
# Набор pdb-фаилов у нас есть после выполнения скрипта dataset_prepare. Нужно подготовить для каждой структуры свой csv с списком мутаций.
# Считаем, что у нас есть dataset_name.csv со следующими столбцами
'''
'PDB' (PDB_id; str), 
'WILD_TYPE' (нативный остаток однобуквенный; str), 
'CHAIN' (цепь, в которой выполняется замена; str),
'POSITION' (номер заменяемого остатка; str), 
'MUTANT' (мутантный остаток однобуквенный; str), 
'DDG' (экспериментально полученное ddG, с ним будем сравнивать и на нем будем учить модель; float)
'PDB_id_chain' (PDB+CHAIN, str)
'''

PDB_dataset = pd.read_csv("/path/to/dataset_name.csv")

# Готовим список фаилов с мутациями для каждого из pdb. Формат строк фаила для каждого из pdb - цепь нативный_остаток+позиция+мутантный_остаток.

for i in PDB_dataset.groupby("PDB"):
    pdb_all_mut = PDB_dataset.groupby("PDB").get_group(i[0])
    pdb_all_mut["mutation_code"] = pdb_all_mut["WILD_TYPE"] + pdb_all_mut["POSITION"] + pdb_all_mut["MUTANT"] #создаем столбец нативный_остаток+позиция+мутантный_остаток
    pdb_all_mut = pdb_all_mut[["CHAIN","mutation_code"]] # оставляем только нужные столбцы
    pdb_all_mut.to_csv(r'/path/to/folder/' + i[0] + ".csv", header=False, index=None, sep=' ', mode='w')

# Отправляем все структуры на веб-инструмент (сейчас - на проверку стабильности). На выходе получаем набор фаилов. 
# Для разных веб-инструментов будут разные команды, в зависимости от html веб-инструмента.

for pdb_id in PDB_dataset['PDB_id_chain'].unique():
    pdb_name = pdb_id[0:5]+"_model_1.pdb" #после dataset_prepare у нас есть набор фаилов такого типа
    mut_name = pdb_id[0:4] + ".csv"
    driver = webdriver.Chrome(ChromeDriverManager().install())
    driver.get('http://biosig.unimelb.edu.au/mcsm/stability') # ссылка на веб-инструмент
    path_pdb = "/path/to/folder/" + pdb_name
    path_mutations = "/path/to/folder/" + mut_name
    driver.find_element(By.XPATH, "/html/body/div[2]/div[2]/div[3]/div/form/input[1]").send_keys(path_pdb) # загружаем pdb фаил
    driver.find_element(By.XPATH, "/html/body/div[2]/div[2]/div[3]/div/form/input[2]").send_keys(path_mutations) # загружаем список мутаций
    driver.find_element(By.XPATH, "/html/body/div[2]/div[2]/div[3]/div/form/button").click() # кликаем на кнопочку "ПОСЧИТАТЬ"
    time.sleep(30) # нужно, чтобы мутации успели посчитаться, можно варьировать
    driver.find_element(By.XPATH, "/html/body/div[2]/div[2]/div/div[2]/div[2]/a").click() # кликаем на кнопочку "сохранить в csv", после чего нас перекинет на html с csv таблицей наших результатов
    data = driver.page_source # сохраняем текущую страницу (с csv)
    data_parsed = BeautifulSoup(data,'html.parser') # парсим html
    # Дальше все зависит от того, как инструмент предоставляет результаты
    # Цель - получить результаты в формате csv со строками, схожими со строками в dataset_name.csv
    
    table = data_parsed.find("pre") 
    result = pd.read_csv(StringIO(table.contents[0]), sep="\t")
    result.to_csv(r'/path/to/folder/' + pdb_id[0:5] + "_model_1_result_web.csv", header=True, index=None, sep=' ', mode='w')
    print("Получены результаты по ", pdb_name) 
    driver.close() # закрываем страницу, чтобы не плодить кучу окон

# Теперь у нас есть набор фаилов, давайте их объединим в один, схожий с dataset_name.csv

PDB_dataset_web_server = pd.DataFrame()
for pdb_id in PDB_dataset['PDB_id_chain'].unique():
    data = pd.read_csv(r'/path/to/folder/' + pdb_id[0:5] + "_model_1_result_web.csv", sep=" ")
    PDB_dataset_web_server = PDB_dataset_web_server.append(data, ignore_index=True)

# Давайте проверим, что наш датасет и полученный датасет после веб-инструмента одинаковой длины

if PDB_dataset_web_server.shape[0] != PDB_dataset.shape[0]:
    print('Что то не так - размер датасета после инструмента не сходится с изначальным датасетом')

# Давайте также проверим, что наш датасет и полученный датасет после веб-инструмента полностью совпадают по своей структуре. Проверьте типы данных для каждого столбца - везде должно быть str.

check = PDB_dataset_web_server[["WILD_TYPE", 'CHAIN', 'POSITION', 'MUTANT']] == PDB_dataset[["WILD_TYPE", 'CHAIN', 'POSITION', 'MUTANT']]
for index, row in check.iterrows():
    if False in row.values:
        print(index)
        print(row)
    
# Если все хорошо и в предыдущем шаге у вас ничего не напечаталось, то можно объединить результаты в один датафрейм

PDB_dataset["DDG_web_tool"] = PDB_dataset_web_server["DDG"]

# И сохраним наш дополненный датасет

PDB_dataset.to_csv(r'/path/to/folder/dataset_name.csv', header=True, mode='w')


import SMILES
import os
from tkinter import filedialog
from tkinter import *
import re
import ast
import string

root = Tk()
root.withdraw()
open_file_selected = filedialog.askopenfile()
smiles_file = open(f'{os.path.abspath(open_file_selected.name).rstrip(".txt")}_poly.txt', 'w+')
num_to_polymers = {1: 'ds', 2: 'db', 3: 'pr', 4: 'pl', 5: 'cou', 6: 'flr', 7: 'tmr', 8: 'c3', 9: 'c5'}
matryca_permutacji = []
for line in open_file_selected:
    do_wymiany = r"\[\'\d-.+\'\]"
    cytat = re.search(do_wymiany, line)
    # matryca_permutacji.append(ast.literal_eval(line[cytat.start():cytat.end()]))
    matryca_permutacji.append(line[cytat.start():cytat.end()])
open_file_selected.close()
nowa_lista = []

for blok in matryca_permutacji:
    nowy_wiersz = ''
    for znak in blok:
        if znak in string.digits:
            nowy_wiersz += num_to_polymers[int(znak)]
        else:
            nowy_wiersz += znak
    nowa_lista.append(ast.literal_eval(nowy_wiersz))

nowa_lista2 = []
for each_blok in range(len(nowa_lista)):
    nowa_lista2.append([])
    for each_polymer in range(len(nowa_lista[each_blok])):
        Smiles_polymer = SMILES.create_polymer(nowa_lista[each_blok][each_polymer])
        nowa_lista2[each_blok].append(Smiles_polymer)

for each in nowa_lista2:
    smiles_file.write(str(each))
smiles_file.close()

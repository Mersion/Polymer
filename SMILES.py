# The program I write here will focus on building a trimer out of any of 8 monomers in smiles.
from pysmiles import read_smiles
import networkx as nx
import matplotlib.pyplot as plt
import string
def create_polymer(sequence):
    n = len(sequence.split('-'))
    seqsplit = sequence.split('-')
    smiles = 'OP(O)(=O)O'
    for i in range(n):
        if i != n - 1:
            if seqsplit[i] in ['c3', 'c5']:
                smiles += polymer_dict[seqsplit[i]] + polymer_dict['h']
                continue
            else:
                smiles += 'CC(%s)C' % (polymer_dict[seqsplit[i]]) + polymer_dict['h']
        else:
            if seqsplit[i] in ['c3', 'c5']:
                smiles += polymer_dict[seqsplit[i]] + 'O'
                continue
            else:
                smiles += 'CC(%s)C' % (polymer_dict[seqsplit[i]]) + 'O'
    return smiles


def reverse(seq):
    switch = 0
    reversedseq = ''
    keeped_seq = ''
    stored_digit = ''
    equal_sign = ''
    for i in range(len(seq), 0, -1):
        if seq[i-1] in [')',']']:
            switch += 1
        if switch > 0:
            keeped_seq += seq[i-1]
        if seq[i-1] in ['(','[']:
            switch -= 1
        if switch == 0 and seq[i-1] in string.printable and seq[i-1] not in ['(','[',']',')']:
            if seq[i-1] in string.digits:
                stored_digit = seq[i-1]
                continue
            if seq[i-1] =='=':
                equal_sign = '='
                continue
            if keeped_seq != '':
                reversedseq += seq[i-1] + stored_digit + f'{equal_sign}' + keeped_seq[::-1]
                keeped_seq = ''
                stored_digit = ''
                equal_sign = ''
            else:
                reversedseq += f'{equal_sign}'+ seq[i-1] + stored_digit
                stored_digit = ''
                equal_sign = ''
    return reversedseq

def reverse_all(dictionary):
    seq = list(polymer_dict.values())
    reversed_seq = []
    for each in seq:
        reversed_seq.append(reverse(each))
    polymer_dictionary = {key: value for key, value in zip(polymer_dict.keys(), reversed_seq)}
    return polymer_dictionary
# nie zmieniamy kolejności wnętrza nawiasów
# gdy napotka numer,nawias przechowuje go do czasu literki, a potem go wrzuca po literce

# Dansyl
ds = 'CN(C)c1cccc2c1cccc2S(=O)(=O)N'
ds_reversed = 'NS(=O)(=O)c2cccc1c2cccc1N(C)C'

# Dabsyl
db = 'CN(C)c1ccc(cc1)N=Nc1ccc(cc1)S(=O)(=O)N'
r_db ='NS(=O)(=O)c(cc1)ccc1N=Nc(cc1)ccc1N(C)C'

# Pyrene
pr = 'C=1CC=2C(C5C1)C4C(C=C5)C=CC(C4=CC2)COC(=O)N'

# Perylene
pl = 'c1ccc2cccc5c2c1c3cccc4c3c5ccc4COC(=O)N'

# Coumarine
cou = 'CCN(C)c1ccc2c(c1)OC(=O)C(C2)C(=O)N'

# Fluoresceina
flr = 'CC(C)(C)C(=O)OC1C=CC=2C4(c3ccc(OC(=O)C(C)(C)C)cc3OC2C=1)OC(=O)C5=C4C=CC=(C5)C(=O)NCCCC'

# Tetrametylorodamina
tmr = 'CN(C)C1C=CC=2C4(c3ccc(N(C)C)cc3OC2C=1)OC(=O)C5=C4C=CC=(C5)OC(=O)N'

# Cyjanina-3
c3 = 'CCCN=3C1=CC=CC=C1C(C)(C)C3C=CC=C1C(C)(C)c2ccccc2N1CCC'

# Cyjanina-5
c5 = 'CCCN=3C1=CC=CC=C1C(C)(C)C3C=CC=CC=C1C(C)(C)c2ccccc2N1CCC'

# Łącznik
alif = 'CCC'

# Kwas Fosforowy H3PO4
h = 'OP(O)(=O)O'

polymer_dict = {'ds': 'NS(=O)(=O)c2cccc1c2cccc1N(C)C', 'db': 'NS(=O)(=O)c(cc1)ccc1N=Nc(cc1)ccc1N(C)C',
                'pr': 'NC(=O)OCC1C=CC2C=CC3C=CCC4=CC=C1C2C34', 'pl': 'NC(=O)OCC1=C2C=CC=C3C4=CC=CC(=C45)C=CC=C5C(=C23)C=C1',
                'cou': 'NC(=O)C(C2)C(=O)Oc(c1)c2ccc1N(CC)CC', 'flr': 'CCCCNC(=O)c1ccc2C(=O)OC3(C4C=CC(=O)=CC=4OC5=CC(=O)=CC=C35)c2c1',
                'tmr': 'CCCCNC(=O)c1ccc2C(=O)OC3(C4C=CC(N(C)C)=CC=4OC5=CC(N(C)C)=CC=C35)c2c1','c3': 'CCCN3C1=CC=CC=C1C(C)(C)C3=C/C=C/C=1C(C)(C)c2ccccc2[N+]1CCC',
                'c5': 'CCCN3C1=CC=CC=C1C(C)(C)C3=C/C=C/C=C/C=1C(C)(C)c2ccccc2[N+]1CCC' , 'alif': 'CCC', 'h': 'OP(O)(=O)O'}



# SMILE = create_polymer('db-ds-pr-pl-cou-c3-flr-c5-tmr')
# print(SMILE)
# print('-'.join(list(polymer_dict.keys())))
# print(reverse(c3))
# print(reverse_all(polymer_dict))
#
#
#
#
#
# SMILE_py = read_smiles(SMILE)
# # atom vector (C only)
# print(SMILE_py.nodes(data='element'))
# # adjacency matrix
# #print(nx.to_numpy_matrix(SMILE_py))
# print(nx.adjacency_matrix(SMILE_py, weight='order').todense())
#
# elements = nx.get_node_attributes(SMILE_py, name = "element")
# nx.draw(SMILE_py, with_labels=True, labels = elements, pos=nx.spring_layout(SMILE_py))
# plt.gca().set_aspect('equal')
# plt.show()
#
# Dansyl = [C12H12ClNO2S, 269.747 ]
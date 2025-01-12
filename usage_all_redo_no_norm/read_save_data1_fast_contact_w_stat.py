import numpy as np
import pandas as pd
import sys, os
from random import shuffle

#from sklearn import metrics

import random

from random import choice




##np.save('aa_vec_dic_300.npy', aalist)

dict_load=np.load('aa_vec_dic_300.npy', allow_pickle=True).item()




def residue_features(resname):
    seq=("ALA","ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",'Unkown')
    return np.array(one_of_k_encoding_unk(resname,seq))



def one_of_k_encoding(x, allowable_set):
    if x not in allowable_set:
        raise Exception("input {0} not in allowable set{1}:".format(x, allowable_set))
    return list(map(lambda s: x == s, allowable_set))

def one_of_k_encoding_unk(x, allowable_set):
    """Maps inputs not in the allowable set to the last element."""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))







aa_codes = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']




aa3to1={
   'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
   'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
   'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
   'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G'
}

aa1to3={value:key  for key,value  in aa3to1.items()}


import sys


frr=open('out_summary.txt','r')
arr_frr=frr.readlines()

all_label=[]
all_env_vec=[]
for line in arr_frr:
    arr_tem=line.split(':::')
    label= arr_tem[0]
    arr_temX=arr_tem[1].split(',')
    envVec= np.zeros(300)

    for residue in arr_temX:
            envVec=envVec+dict_load[residue.strip()]
    all_label.append(label)
    all_env_vec.append(envVec)


np.save('all_label.npy',all_label)
np.save('all_env_vec.npy',all_env_vec)




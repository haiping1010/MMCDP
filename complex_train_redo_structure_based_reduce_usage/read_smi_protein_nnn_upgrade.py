import os
import sys
import numpy as np
import glob
from torch_geometric.data import InMemoryDataset, DataLoader
from torch_geometric import data as DATA
import torch

from rdkit import Chem
from rdkit.Chem import MolFromSmiles
import networkx as nx
import pandas as pd
import re, os
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
from MDAnalysis.analysis import distances

import random
from torch import nn





#idx=sys.argv[1]


def atom_features(atom):
    return np.array(one_of_k_encoding_unk(atom.GetSymbol(),['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na','Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb','Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'H','Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr','Cr', 'Pt', 'Hg', 'Pb', 'Unknown']) +
                    one_of_k_encoding(atom.GetDegree(), [0, 1, 2, 3, 4, 5, 6,7,8,9,10]) +
                    one_of_k_encoding_unk(atom.GetTotalNumHs(), [0, 1, 2, 3, 4, 5, 6,7,8,9,10]) +
                    one_of_k_encoding_unk(atom.GetImplicitValence(), [0, 1, 2, 3, 4, 5, 6,7,8,9,10]) +
                    [atom.GetIsAromatic()])




def one_of_k_encoding(x, allowable_set):
    print (x,allowable_set)
    if x not in allowable_set:
        raise Exception("input {0} not in allowable set{1}:".format(x, allowable_set))
    return list(map(lambda s: x == s, allowable_set))

def one_of_k_encoding_unk(x, allowable_set):
    """Maps inputs not in the allowable set to the last element."""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))




def calc_atom_features(atom, explicit_H=False):
    """
    atom: rdkit.Chem.rdchem.Atom
    explicit_H: whether to use explicit H
    use_chirality: whether to use chirality
    """
    results = one_of_k_encoding_unk(
      atom.GetSymbol(),
      [
       'C', 'N', 'O', 'S', 'F', 'P', 'Cl',
                'Br', 'I', 'B', 'Si', 'Fe', 'Zn',
                'Cu', 'Mn', 'Mo', 'other'
      ]) + one_of_k_encoding(atom.GetDegree(),
                             [0, 1, 2, 3, 4, 5, 6,7,8,9,10]) + \
              [atom.GetFormalCharge(), atom.GetNumRadicalElectrons()] + \
              one_of_k_encoding_unk(atom.GetHybridization(), [
                Chem.rdchem.HybridizationType.SP, Chem.rdchem.HybridizationType.SP2,
                Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP3D,
                Chem.rdchem.HybridizationType.SP3D2,'other']) + [atom.GetIsAromatic()]
                # [atom.GetIsAromatic()] # set all aromaticity feature blank.
    # In case of explicit hydrogen(QM8, QM9), avoid calling `GetTotalNumHs`
    if not explicit_H:
        results = results + one_of_k_encoding_unk(atom.GetTotalNumHs(),
                                                  [0, 1, 2, 3, 4])
    return np.array(results)


def calc_bond_features(bond, use_chirality=True):
    """
    bond: rdkit.Chem.rdchem.Bond
    use_chirality: whether to use chirality
    """
    bt = bond.GetBondType()
    bond_feats = [
        bt == Chem.rdchem.BondType.SINGLE, bt == Chem.rdchem.BondType.DOUBLE,
        bt == Chem.rdchem.BondType.TRIPLE, bt == Chem.rdchem.BondType.AROMATIC,
        bond.GetIsConjugated(),
        bond.IsInRing()
    ]
    if use_chirality:
        bond_feats = bond_feats + one_of_k_encoding_unk(
            str(bond.GetStereo()),
            ["STEREONONE", "STEREOANY", "STEREOZ", "STEREOE"])
    return np.array(bond_feats+[False] * 10).astype(int)




def load_mol(molpath, explicit_H=False, use_chirality=True):
        # load mol
        if re.search(r'.pdb$', molpath):
                mol = Chem.MolFromPDBFile(molpath, removeHs=not explicit_H)
        elif re.search(r'.mol2$', molpath):
                mol = Chem.MolFromMol2File(molpath, removeHs=not explicit_H)
        elif re.search(r'.sdf$', molpath):
                mol = Chem.MolFromMolFile(molpath, removeHs=not explicit_H)
        else:
                raise IOError("only the molecule files with .pdb|.sdf|.mol2 are supported!")

        if use_chirality:
                Chem.AssignStereochemistryFrom3D(mol)
        return mol



import pickle
aa_dict=np.load('aa_vec_dic.npy', allow_pickle=True).item()
cord = [None] * 3

with open('all_data/dict_poc_pos.pickle', 'rb') as f:
    dict_poc_loaded = pickle.load(f)

print (dict_poc_loaded.keys())

def  pdb_graph(antib_file,antig_file):
  uniq=[]
  Pposition={}
  ResinameP={}
  Interface=[]
  residuePair=[]
  PL_pair=[]
  #arr_all=dict_poc_loaded[pdbfile]
  #residuePair,uniq,Pposition=arr_all[0],arr_all[1],arr_all[2]
  residuePair_P,uniq_P,Pposition,ResinameP=dict_poc_loaded[antib_file]
  residuePair_E,uniq_E,Eposition,ResinameE=dict_poc_loaded[antig_file]

  for key1, value1 in Pposition.items():
      for key2, value2 in Eposition.items():
            a = np.array(value1)
            a1 = a.astype(np.float)
            b = np.array(value2)
            b1 = b.astype(np.float)
            xx=np.subtract(a1,b1)
            tem=np.square(xx)
            tem1=np.sum(tem)
            out=np.sqrt(tem1)
            if out<9:
                PL_pair.append([ResinameP[key1],ResinameE[key2],out])
                uniq.append(ResinameP[key1])
                #uniq.append('mol'+index)
                #Interface.append(a1)

  uniq_n=list(set(uniq_P))
  my_dict = {}
  for index, item in enumerate(uniq_n):
        my_dict[item] = index

  uniq_n_E=list(set(uniq_E))
  my_dict_E = {}
  for index, item in enumerate(uniq_n_E):
        my_dict_E[item] = index


  edges_p=[]
  features=[]
  edge_attrs=[]
  for i in residuePair_P:
     edges_p.append([my_dict[i[0]], my_dict[i[1]]])
     code_d=int(i[2])
     edge_attr=one_of_k_encoding(code_d,[0, 1, 2, 3, 4, 5, 6,7,8,9])
     edge_attr=edge_attr
     edge_attrs.append(edge_attr)

  for i in PL_pair:
     print(my_dict_E.keys())
     print (i[0],i[1])
     if i[0] in my_dict.keys()  and i[1] in my_dict_E.keys():
        edges_p.append([my_dict[i[0]], int(my_dict_E[i[1]]+len(uniq_n))])
     
        code_d=int(i[2])
        edge_attr=one_of_k_encoding(code_d,[0, 1, 2, 3, 4, 5, 6,7,8,9])
        edge_attr=edge_attr
        print (len(edge_attr))
        edge_attrs.append(edge_attr)

  for i in residuePair_E:
     edges_p.append([int(my_dict_E[i[0]]+len(uniq_n)), int(my_dict_E[i[1]]+len(uniq_n))])
     code_d=int(i[2])
     edge_attr=one_of_k_encoding(code_d,[0, 1, 2, 3, 4, 5, 6,7,8,9])
     edge_attr=edge_attr
     edge_attrs.append(edge_attr)

  print (len(edge_attrs))
  
  for index, item in enumerate(uniq_n):
        #print (item)
        feature = aa_dict[item[0:3]]
        feature = feature.tolist() + [0] * (60 - len(feature))
        #feature.append(1)
        feature=np.array(feature).astype(int)
        features.append( feature )

        print (len(feature))

  for index, item in enumerate(uniq_n_E):
        #print (item)
        feature = aa_dict[item[0:3]]
        feature = [0] * (60 - len(feature)) + feature.tolist() 
        #feature.append(1)
        feature=np.array(feature).astype(int)
        features.append( feature )

        print (len(feature))

  c_size=len(uniq_n)+len(uniq_n_E)
  
  #print (c_size)
  print (len(edge_attrs),'gggggggggggggggg')
  print (len(edges_p))
  return  c_size,features,edges_p, edge_attrs

pocket_graph = {}
smile_graph={}
import glob

frr=open('all_data/'+'temT.txt','r')
arr_frr=frr.readlines()
pro_lig={}
arr_name=[]


import json
# Load dict_label from JSON file
#with open('all_data/dict_aff.json', 'r') as f:
#    dict_label_loaded = json.load(f)

train_Y=[]
for name in arr_frr:
    #mol = Chem.MolFromPDBFile('all_data/'+name.strip()+'_ligand_n.pdb', removeHs=not False)
    #mol=load_mol('all_data/'+name[0:4]+'_ligand_fixed.pdb')
    #if  mol is not None:
    arr_name.append(name.strip())
    #train_Y.append(round(float(dict_label_loaded[name.strip()]),2))
    



#df = pd.read_csv('data/' + dataset + '_train.csv')
#train_drugs, train_prots,  train_Y = list(df['compound_iso_smiles']),list(df['target_sequence']),list(df['affinity'])
#XT = [seq_cat(t) for t in train_prots]
#train_drugs, train_prots,  train_Y = np.asarray(train_drugs), np.asarray(XT), np.asarray(train_Y)

#train_data = TestbedDataset(root='data1', dataset='pocket_train', xd=arr_name, smile_graph=pocket_graph)

class TestbedDataset(InMemoryDataset):
    def __init__(self, root='/tmp', dataset='davis', 
                 xd=None, pocket_graph=None, y=None, transform=None,
                 pre_transform=None,smile_graph=None):

        #root is required for save preprocessed data, default is '/tmp'
        super(TestbedDataset, self).__init__(root, transform, pre_transform)
        # benchmark dataset, default = 'davis'
        self.dataset = dataset
        #print (self.processed_paths[0]+'xxxxxxxxxxx')
        self.process(xd,pocket_graph,y,smile_graph)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        pass
        #return ['some_file_1', 'some_file_2', ...]

    @property
    def processed_file_names(self):
        return [self.dataset + '.pt']

    def download(self):
        # Download to `self.raw_dir`.
        pass

    def _download(self):
        pass

    def _process(self):
        if not os.path.exists(self.processed_dir):
            os.makedirs(self.processed_dir)

    # Customize the process method to fit the task of drug-target affinity prediction
    # Inputs:
    # XD - list of SMILES, XT: list of encoded target (categorical or one-hot),
    # Y: list of labels (i.e. affinity)
    # Return: PyTorch-Geometric format processed data
    def process(self, xd, pocket_graph, y, smile_graph):
        
        data_list = []
        data_len = len(xd)
        for i in range(data_len):
            print('Converting SMILES to graph: {}/{}'.format(i+1, data_len))
            filename = xd[i]
            #print (filename)
            labels=y[i]
            # convert SMILES to molecular representation using rdkit
            # make the graph ready for PyTorch Geometrics GCN algorithms:
            #c_size, features, edge_index =  smile_graph[str(filename[0:4])]
            #print (smile_graph[str(filename[0:4])])
            print (filename)
            tem_name=filename.split()
            c_size, features, edge_index,edge_attrs =  pdb_graph(tem_name[0].replace('.pdb','_poc.pdb'), tem_name[1].replace('.pdb','_poc.pdb'))
            #print (smile_graph[str(filename[0:4])])
            GCNData = DATA.Data(x=torch.Tensor(features),
                                edge_index=torch.LongTensor(edge_index).transpose(1, 0),
                                edge_attr=torch.Tensor(edge_attrs),
                                y=torch.FloatTensor([labels]))
            GCNData.__setitem__('c_size', torch.LongTensor([c_size]))

            #####
          
            #print (pocket_graph[str(filename[0:4])])
            GCNData.name = tem_name[1].replace('_poc.pdb','')
            
            # append graph, label and target sequence to data list
            data_list.append(GCNData)


        print('Graph construction done. Saving to file.')
        data, slices = self.collate(data_list)
        # save preprocessed data:
        torch.save((data, slices), self.processed_paths[0])


train_Y=np.ones(len(arr_name))
#train_Y=np.array(train_Y)

train_data = TestbedDataset(root='data1', dataset='L_P_train', xd=arr_name, pocket_graph=pocket_graph, y=train_Y, smile_graph=smile_graph)




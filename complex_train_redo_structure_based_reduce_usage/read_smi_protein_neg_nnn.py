import os
import sys
import numpy as np
from sklearn.cluster import KMeans
import numpy as np
from sklearn.cluster import KMeans
import glob
from torch_geometric.data import InMemoryDataset, DataLoader
from torch_geometric import data as DATA
import torch

from rdkit import Chem
from rdkit.Chem import MolFromSmiles
import networkx as nx

idx=sys.argv[1]

aa_dict=np.load('aa_vec_dic.npy', allow_pickle=True).item()

cord = [None] * 3

Pposition={}
ResinameP={}
Interface=[]
residuePair=[]


uniq=[]


def  pdb_graph(pdbfile):
  uniq=[]
  Pposition={}
  ResinameP={}
  Interface=[]
  residuePair=[]
  for line in open(pdbfile):
     tem_B=' '
     if len(line)>16:
        tem_B=line[16]
        line=line[:16]+' '+line[17:]
     #print(line)
     list_n = line.split()
     id = list_n[0]
     if id == 'ATOM' and tem_B !='B' and line.find(" HOH ") == -1:
        type = list_n[2]
        #print (line)
        if type == 'CA' and list_n[3]!= 'UNK':
            residue = list_n[3]
            atomname=list_n[2]
            type_of_chain = line[21:22]
            tem1=line[22:26].replace("A", "")
            tem2=tem1.replace("B", "")
            tem2=tem2.replace("C", "")

            #tem2=filter(str.isdigit, list_n[5])
            #atom_count = tem2+line[21:22]
            atom_count = line[4:11]+line[21:22]
            cord[0]=line[30:38]
            cord[1]=line[38:46]
            cord[2]=line[46:54]
            position = cord[0:3]
            Pposition[atom_count]=position
            ResinameP[atom_count]=line[17:26]
            #print atom_count,hash[residue[0:3]+atomname]

  for key1, value1 in Pposition.items():
     for key2, value2 in Pposition.items():
         if key2>key1:
            a = np.array(value1)
            a1 = a.astype(np.float)
            b = np.array(value2)
            b1 = b.astype(np.float)
            xx=np.subtract(a1,b1)
            tem=np.square(xx)
            tem1=np.sum(tem)
            out=np.sqrt(tem1)
            if out<5 :
                residuePair.append([ResinameP[key1],ResinameP[key2]])
                uniq.append(ResinameP[key1])
                uniq.append(ResinameP[key2])
                Interface.append(a1)
  uniq_n=list(set(uniq))
  my_dict = {}
  for index, item in enumerate(uniq_n):
        my_dict[item] = index

  edges_p=[]
  features=[]
  for i in residuePair:
     #print ( my_dict[i[0]], my_dict[i[1]])
     edges_p.append([my_dict[i[0]], my_dict[i[1]]])
  for index, item in enumerate(uniq_n):
        #print (item)
        feature = aa_dict[item[0:3]]
        features.append( feature )
        #features.append( feature / sum(feature) )
  c_size=len(uniq_n)
  #print (c_size)
  return  c_size,features,edges_p

pocket1_graph = {}
pocket2_graph={}
import glob

frr=open('tem_'+str(idx)+'.txt','r')

##frr=open('list_final_nn.txt','r')
arr_frr=frr.readlines()
pro_lig={}
arr_name=[]



for name in arr_frr:
    arr_name.append(name.strip())


#df = pd.read_csv('data/' + dataset + '_train.csv')
#train_drugs, train_prots,  train_Y = list(df['compound_iso_smiles']),list(df['target_sequence']),list(df['affinity'])
#XT = [seq_cat(t) for t in train_prots]
#train_drugs, train_prots,  train_Y = np.asarray(train_drugs), np.asarray(XT), np.asarray(train_Y)

#train_data = TestbedDataset(root='data1', dataset='pocket_train', xd=arr_name, smile_graph=pocket_graph)

class TestbedDataset(InMemoryDataset):
    def __init__(self, root='/tmp', dataset='davis', 
                 xd=None, pocket1_graph=None, y=None, transform=None,
                 pre_transform=None,pocket2_graph=None):

        #root is required for save preprocessed data, default is '/tmp'
        super(TestbedDataset, self).__init__(root, transform, pre_transform)
        # benchmark dataset, default = 'davis'
        self.dataset = dataset
        #print (self.processed_paths[0]+'xxxxxxxxxxx')
        self.process(xd,pocket1_graph,y,pocket2_graph)
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
    def process(self, xd, pocket1_graph, y, pocket2_graph):
        
        data_list = []
        data_len = len(xd)
        for i in range(data_len):
            print('Converting SMILES to graph: {}/{}'.format(i+1, data_len))
            filename = xd[i]
            #print (filename)
            labels=y[i]
            print (filename)
            arr_f=filename.split()
            # convert SMILES to molecular representation using rdkit
            # make the graph ready for PyTorch Geometrics GCN algorithms:
            #c_size, features, edge_index =  smile_graph[str(filename[0:4])]
            #fr=open(str(arr_f[0])+'_poc.pdb','r')
            #print (linearr[0])
            c_size, features, edge_index =  pdb_graph(arr_f[0].replace('.pdb','_poc.pdb'))
            #print (smile_graph[str(filename[0:4])])
            GCNData = DATA.Data(x=torch.Tensor(features),
                                edge_index=torch.LongTensor(edge_index).transpose(1, 0),
                                y=torch.FloatTensor([labels]))
            GCNData.__setitem__('c_size', torch.LongTensor([c_size]))
          
            c_size1, features1, edge_index1 =  pdb_graph(arr_f[1].replace('.pdb','_poc.pdb'))
            #print (pocket_graph[str(filename[0:4])])
            GCNData.name = arr_f[0].replace('antigen.pdb','')+"-"+arr_f[1].replace('antibody.pdb','')
            GCNData.target = DATA.Data(x=torch.Tensor(features1),
                                edge_index=torch.LongTensor(edge_index1).transpose(1, 0))
            
            # append graph, label and target sequence to data list
            data_list.append(GCNData)


        print('Graph construction done. Saving to file.')
        data, slices = self.collate(data_list)
        # save preprocessed data:
        torch.save((data, slices), self.processed_paths[0])


train_Y=np.zeros(len(arr_name))
train_data = TestbedDataset(root='data1', dataset='L_P_train_neg_'+idx, xd=arr_name, pocket1_graph=pocket1_graph, y=train_Y, pocket2_graph=pocket2_graph)




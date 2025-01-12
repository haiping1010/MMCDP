# -*- coding: utf-8 -*-
import sys

filebase=sys.argv[1]


def readpdb(filename):
   Pposition={}
   ResinameP={}
   Atomline=[]
   for line in open(filename):
   #for line in open(filebase+'_antibody.pdb'):
      tem_B=' '
      if len(line)>16:
         tem_B=line[16]
         line=line[:16]+' '+line[17:]
      #print(line)
      list_n = line.split()
      id = list_n[0]
      if id == 'ATOM' and tem_B !='B':
         type = list_n[2]
         #if type == 'CA' and list[3]!= 'UNK':
         if   list_n[3]!= 'UNK' :
             residue = list_n[3]
             type_of_chain = line[21:22]
             tem1=line[22:26].replace("A", "")
             tem2=tem1.replace("B", "")
             tem2=tem2.replace("C", "")

             #tem2=filter(str.isdigit, list[5])
             atom_count = line[4:11]+tem2+line[21:22]
             list_n[6]=line[30:38]
             list_n[7]=line[38:46]
             list_n[8]=line[46:54]
             position = list_n[6:9]
             Pposition[atom_count]=position
             list_n[4]=line[21:22]
             list_n[5]=line[22:27].replace(' ','')
             ResinameP[atom_count]=residue+list_n[5]+list_n[4]
             resindex=residue+list_n[5]
             #Atomline.append(line)
             Atomline.append(line)
             Pposition[atom_count]=position
   return   Pposition,ResinameP,Atomline

def save_pdb(pdb_all_line,name):
     fw=open(name,"w")
     for line in pdb_all_line:
         fw.write(line)
     fw.close()

HLposition, ResinameHL, Atomline=readpdb(filebase+'_antibody.pdb')
Eposition, ResinameE,AtomlineE=readpdb(filebase+'_antigen.pdb')





import itertools
from collections import OrderedDict


import pickle
from collections import OrderedDict

import pickle
import glob
import numpy as np
possible_mutations = []


for filename in glob.glob('identify_point_*_possible_mutations.pkl'):
#for filename in glob.glob('identify_point*.pkl'):
    with open(filename, 'rb') as f:
         data = pickle.load(f)
         #print (data.keys())
         possible_mutations= possible_mutations + list(data.keys())

print(possible_mutations)



all_pair=[]
tem_neighbor=[]
HL_res=[]
for key1, value1 in HLposition.items():
     #print key1[7:].replace(' ','') , possible_mutations
     if key1[7:].replace(' ','') in possible_mutations:
        for key2, value2 in Eposition.items():
            a = np.array(value1)
            a1 = a.astype(np.float)
            b = np.array(value2)
            b1 = b.astype(np.float)
            xx=np.subtract(a1,b1)
            tem=np.square(xx)
            tem1=np.sum(tem)
            #out=np.sqrt(tem1)
            #print out
            if tem1<36 :
                tem_neighbor.append([ResinameHL[key1], ResinameE[key2]])
                HL_res.append(ResinameHL[key1])
                #print out
   #print  (ResinameHL[key1],tem_neighbor)

   #fo = open(filex+'_'+filey+"_summary.txt", "w")

#tem_neighbor=list(set(tem_neighbor))
unique_tuples = set(tuple(item) for item in tem_neighbor)
tem_neighbor = [list(item) for item in unique_tuples]
HL_res=list(set(HL_res))
print (HL_res)
for  linex in HL_res:
    print(linex)
    index_nn=0
    for liney in tem_neighbor:
        #print (linex, liney[0])
        if  linex==liney[0]:
            index_nn=index_nn+1

            #print (','.join(tem_neighbor))
            if index_nn==1:
                str1=linex+':::'+liney[1][0:3]
            else:
                str1=str1+','+liney[1][0:3]
    all_pair.append(str1+"\n")

fo = open("out_summary.txt", "w")
for name in all_pair:
    fo.write(name)

fo.close()











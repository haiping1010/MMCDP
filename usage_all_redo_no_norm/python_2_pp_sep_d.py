import sys
import numpy as np
from sklearn.cluster import KMeans
import numpy as np
from sklearn.cluster import KMeans
def one_hot(i):
    m = np.zeros(20, 'uint8')
    m[i] = 1
    return m

def MapResToOnehot(residue1,residue2):
   code={}
   seq=("ALA","ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
   index=0
   for name in seq:
              #print (name)
              code[name]=one_hot(index)
              index=index+1
   return (np.concatenate((code[residue1],code[residue2])))
#print (MapResToOnehot("ARG","VAL"))


      
HLposition={}
Eposition={}
ResinameHL={}
ResinameE={}
residuePair=[]
#if len(sys.argv) <1 :
#   print("python python2.py xxx.pdb")
#file=sys.argv[1]
#filebase=file
import glob

def  read_pdb(pdbname):
  Atomline=[]
  HLposition={}
  ResinameHL={}

  for line in open(pdbname):
     #Atomline.append(line)
     tem_B=' '
     if len(line)>16:
        tem_B=line[16]
        line=line[:16]+' '+line[17:]
     #print(line)
     list = line.split()
     id = list[0]
     if id == 'ATOM' and tem_B !='B':
        Atomline.append(line)
        type = list[2]
        if type == 'CA'  and list[3]!= 'UNK'  :
            residue = list[3]
            type_of_chain = line[21:22]
            ##tem1=list[5].replace("A", "")
            tem1=line[22:26].replace("A", "")
            tem2=tem1.replace("B", "")
            tem2=tem2.replace("C", "")

            #tem2=filter(str.isdigit, list[5])
            atom_count = tem2+line[21:22]
            list[6]=line[30:38]
            list[7]=line[38:46]
            list[8]=line[46:54]
            position = list[6:9]
            HLposition[atom_count]=position
            list[4]=line[21:22]
            list[5]=line[22:26].replace(' ','')
            ResinameHL[atom_count]=residue+list[5]+list[4]
  return   HLposition, ResinameHL, Atomline



def  caculate_pair(HLposition,Eposition,ResinameHL,ResinameE):
    residuePair=[]
    for key1, value1 in HLposition.items():
       #print (ResinameHL[key1], 'corresponds to', value1)
       for key2, value2 in Eposition.items():
            #print key2,value2
            #print (ResinameE[key2], 'corresponds to', value2)
            ##distance=pow(value1[0]-value2[0])
            a = np.array(value1)
            a1 = a.astype(np.float)
            b = np.array(value2)
            b1 = b.astype(np.float)
            xx=np.subtract(a1,b1) 
            tem=np.square(xx)
            tem1=np.sum(tem)
            out=np.sqrt(tem1)
            #print a
            if out<15 :
                #print (ResinameHL[key1],ResinameHL[key2])
                #print (a,b,out)
                #print (a1)              
                residuePair.append([ResinameHL[key1],ResinameE[key2]])
                #print (antiface)              
    return   residuePair            


import sys

input_f=sys.argv[1]
fr_inf=open(input_f,"r")
arr_inf=fr_inf.readlines()

fw=open(input_f+"_contact.txt","w")
for line in arr_inf:
    tem_arr=line.split()
    antig=tem_arr[0]
    anti=tem_arr[1]
    HLposition, ResinameHL, AtomlineHL = read_pdb(anti)
    Eposition, ResinameE,AtomlineE= read_pdb(antig)

    residuePair=caculate_pair(HLposition,Eposition,ResinameHL,ResinameE)
    
    foo = open(anti.replace(".pdb","_poc.pdb"), "w")

    length=len(residuePair)
    pdbid=antig+"-"+anti
    fw.write(pdbid+","+","+ str(length))
    fw.write('\n')
    for  value1 in AtomlineHL:
        for value2 in residuePair:
           list2 = value1.split()
           #print (value1)
           list2[4]=value1[21:22]
           list2[5]=value1[22:26].replace(' ','')
           resnameId=list2[3]+list2[5]+list2[4]
           #print value2, resnameId
           if  value2[0] == resnameId:
              #print (value2[0])
              foo.write(value1 )
              break
    foo.close()
    foo = open(antig.replace(".pdb","_poc.pdb"), "w")

    for  value1 in AtomlineE:
        for value2 in residuePair:
           list2 = value1.split()
           list2[4]=value1[21:22]
           list2[5]=value1[22:26].replace(' ','')
           resnameId=list2[3]+list2[5]+list2[4]
           #print value2, resnameId
           #print ( value2[1])
           if  value2[1] == resnameId:
              #print (value2[1])
              foo.write(value1 )
              break
    foo.close()



  
fw.close()

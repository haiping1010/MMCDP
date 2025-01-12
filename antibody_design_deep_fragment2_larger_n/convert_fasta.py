import sys
input_f=sys.argv[1]

fr=open(input_f,'r')
fw=open(input_f.replace('.txt','.fasta'),'w')
arr_tem=fr.readlines()
for  idx in range(len(arr_tem)):
     fw.write('>id_'+str(idx).ljust(4,'0')+"\n")
     fw.write(arr_tem[idx].strip()+"\n")
fw.close()

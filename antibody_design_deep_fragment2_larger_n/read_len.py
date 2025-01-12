
import os
import sys

input_f=sys.argv[1]

fr=open(input_f,'r')
fw=open(input_f.replace('fasta','fastax'),'w')
all_arr=fr.readlines()
for name in all_arr:
   if name.startswith('>'):
        oldname=name
   if not name.startswith('>'):
       if len(name.strip()) >3:
          fw.write(oldname)
          fw.write(name)
       else:
          print (name)
          
        
fr.close()
fw.close()

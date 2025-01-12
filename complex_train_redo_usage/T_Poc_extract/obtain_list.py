import shutil, random, os
#filenames = random.sample(os.listdir(dirpath), 4000)
import glob
##filenames = random.sample(os.listdir(dirpath), 3)
#filenames = random.sample(aa, 3)

bb=glob.glob("new_antibody_*_f.pdb")



fw=open("temT.txt", "w")
print (fw.name)
pname='xxxx_antigen.pdb'
for lname in bb:
            fw.write(pname+" "+ lname+"\n")

lname='xxxx_antibody.pdb'
fw.write(pname+" "+ lname+"\n")
fw.close()

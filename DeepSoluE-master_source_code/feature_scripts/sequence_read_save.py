#!/usr/bin/env python
#_*_coding:utf-8_*_

import re, os, sys

def read_protein_sequences(file):
    if os.path.exists("./sequence/"+file) == False:
        print('Error: "' + file + '" does not exist.')
        sys.exit(1)

    with open("./sequence/"+file) as f:
        records = f.read()

    if re.search('>', records) == None:
        print('The input file seems not in fasta format.')
        sys.exit(1)

    records = records.split('>')[1:]
    fasta_sequences = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
        fasta_sequences.append([name, sequence])
    return fasta_sequences
    
    
def save_to_csv(encodings, file):
    with open("./features/"+file, 'w') as f:
        for line in encodings[1:]:
            f.write(str(line[0]))
            for i in range(1,len(line)):
                f.write(',%s' % line[i])
            f.write('\n')

def file_remove():
    import os,shutil
    dir_list1=os.listdir("./")
    for x in dir_list1:
        if x.split("_")[0]=="TMHMM":
             dir_list11=os.listdir("./"+str(x))
             for x1 in dir_list11:
                     os.remove("./"+str(x)+"/"+str(x1))
             shutil.rmtree("./"+str(x))

    dir_list2=os.listdir("./features/")
    for x in dir_list2:
        if x.split(".")[-1]=="csv":
            os.remove("./features/"+str(x))
    
    dir_list3=os.listdir("./features/biofea/")
    for x in dir_list3:
        #if x.split(".")[-1]=="csv":
            os.remove("./features/biofea/"+str(x))

    dir_list4=os.listdir("./sequence/")
    for x in dir_list4:
          if x.split(".")[-1]=="txt":
              os.remove("./sequence/"+str(x))
 
    # dir_list5=os.listdir("./results/")
    # for x in dir_list5:
        # os.remove("./results/"+str(x))
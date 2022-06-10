# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 22:40:49 2020

@author: Administer
"""

import numpy as np
import pandas as  pd
from Bio import SeqIO 




def kmer_seq(sequence, k):
    list_seq=list(sequence)
    return (list("".join(list_seq[i:i+k]) for i in range(len(sequence)-k+1)))


def w2v_kmer_corpus(input_file):
    with open('./sequence/w2v_kmer.txt',"w") as fil1:
            for seq_record in SeqIO.parse('./sequence/'+str(input_file), "fasta"):
                seq_id=seq_record.id
                seq=seq_record.seq
                sub_list=kmer_seq(seq,3)
                fil1.write(str(seq_id)+",")
                fil1.writelines(str(x)+"," for x in sub_list)
                fil1.write("\n")

def word_vector(sequences,model):
    word_features=np.zeros((len(sequences), 100))
    for i in range(len(sequences)):
        array=np.array([model.wv[w] for w in sequences[i] if w in model.wv ])
        idx = np.argwhere(np.all(array[..., :] == 0, axis=1))
        array = np.delete(array, idx, axis=0)
        S= pd.Series(array.mean(axis=0))
        word_features[i]=S.values
    return word_features


def w2v_features():
    import pandas as pd
    import numpy as np
    from gensim import utils
    from gensim.models.word2vec import Word2Vec
    
    with utils.smart_open('./sequence/w2v_kmer.txt','r',encoding='utf-8-sig',) as infile:
        its_list=list(infile)
    list_seq=[] # 保存 kmer序列
    for x in range(len(its_list)): 
        list_seq.append(its_list[x].split(",")[1:])
    
    model= Word2Vec.load("./model/w2v/training_k3w2_shffule.model")       
    X=word_vector(np.array(list_seq),model)
    pd.DataFrame(X).to_csv("features/w2v_feature.csv",sep=",")   
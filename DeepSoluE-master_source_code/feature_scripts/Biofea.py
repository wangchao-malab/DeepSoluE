# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 15:55:39 2021

@author: Administer
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import argparse
from Bio.SeqUtils import ProtParam
from Bio import SeqIO
import numpy as np
import csv, os
import pandas as pd
# Kyte & Doolittle index of hydrophobicity
kd = {'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
      'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
      'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
      'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2}
IUPACProtein_letters='ACDEFGHIKLMNPQRSTVWY'
def check_seq(seq):
    for aa in seq:
        if aa not in IUPACProtein_letters:
            return False
    return True

def fracnumcharge(aa_freq):
    return aa_freq['R'] + aa_freq['K'] + aa_freq['D'] + aa_freq['E']

def thermostability(aa_freq):
    return aa_freq['I'] + aa_freq['V'] + aa_freq['Y'] + aa_freq['W'] + \
           aa_freq['R'] + aa_freq['E'] + aa_freq['L']

def kr_ratio(aa_freq):
    if aa_freq['R'] == 0:
        return np.nan
    else:
        return aa_freq['K']/aa_freq['R']

def de_mul(aa_freq):
    return aa_freq['D']*aa_freq['E']

features = {"sid":[],
                "fracnumcharge": [], "kr_ratio":[],
                "aa_helix": [], "aa_sheet":[], "aa_turn": [],
                "molecular_weight": [], "length": [],
                'avg_molecular_weight': [], "aromaticity": [],
                "instability_index": [], "flexibility": [],
                "gravy": [], "isoelectric_point": [],}


# biopython_features num=12
def biopython_12(input_file):
    with open('./features/biofea/biofea.csv', "w",newline="") as csv_out:
        csv_wr = csv.DictWriter(csv_out, fieldnames=features)
        csv_wr.writeheader()
        for seq in SeqIO.parse("./sequence/"+str(input_file), "fasta"):
            row = dict.fromkeys(features)
            row['sid'] = seq.id
            analysis = ProtParam.ProteinAnalysis(str(seq.seq))
            aa_freq = analysis.get_amino_acids_percent()# Calculate the amino acid content in percentages.
            
            row['fracnumcharge'] = fracnumcharge(aa_freq)
            row['length'] = analysis.length
            row['kr_ratio'] = kr_ratio(aa_freq)
           
            row['molecular_weight'] = analysis.molecular_weight()
            
            row['avg_molecular_weight'] = row['molecular_weight']/row['length']
            row['aromaticity'] = analysis.aromaticity()
            row['instability_index'] = analysis.instability_index()
            row['flexibility'] = np.mean(analysis.flexibility())
            row['gravy'] = analysis.gravy()
            row['isoelectric_point'] = analysis.isoelectric_point()
            h, s, t = analysis.secondary_structure_fraction()
            row["aa_helix"] = h
            row["aa_sheet"] = s
            row["aa_turn"] = t
            csv_wr.writerow(row)




def scale(x):
    import pandas as pd
    import numpy as np
    df1=pd.read_csv("./data/"+str(x),index_col=None,header=None).to_numpy()
    scale_ori={}
    for i in range(len(df1)):
        s=df1[i][0].split(": ")
        scale_ori[str(s[0])]=s[1].strip("  ")
    scale_aa = ['A', 'R', 'N', 'D', 'C','Q', 'E', 'G', 'H', 'I','L', 'K', 'M', 'F', 'P','S', 'T', 'W', 'Y', 'V']
    scale_new={}

    from Bio.SeqUtils import seq3
    for x in scale_aa:
        scale_new[x]=np.float(scale_ori[seq3(x)])
    return scale_new
    



def physico_chemical(fa_path, f_csv,scale_name,w,xx):
    aa_scale=scale(scale_name)
    #print(aa_scale)
    features = {"sid":[],
                }
    for window in [xx]:
        features[str(window)]=[]
        
    with open("./features/biofea/"+f_csv, "w",newline="") as csv_out:
        csv_wr = csv.DictWriter(csv_out, fieldnames=features)
        csv_wr.writeheader()
        for seq in SeqIO.parse(fa_path, "fasta"):
            if not check_seq(seq.seq):
                continue
            row = dict.fromkeys(features)
            row['sid'] = seq.id
            analysis = ProtParam.ProteinAnalysis(str(seq.seq))
                
            f_scale = np.mean(analysis.protein_scale(aa_scale,w))
            row[str(xx)] = f_scale
            csv_wr.writerow(row)

 
def scale_main(input_file):
    i=0
    scalew=[3,7,5,]
    for x in ['14_hc.txt','38_hj.txt','43_hh.txt']:
        #print(x)
        xx=x.split("_")[0]
        w=scalew[i]
        physico_chemical("./sequence/"+str(input_file), xx+".csv",x,w,xx)
        i=i+1


# usearch feature extract
def process_blast6(blast6, chosen):
    names = ("query", "target", "identity", "alignment_length",
             "number_of_mismatches", "number_of_gap_opens",
             "start_position_in_query",  "end_position_in_query",
             "start_position_in_target", "end_position_in_target",
             "e-value", "bit_score")
    b6 = pd.read_csv("./features/biofea/"+blast6, names=names, sep="\t")
    if chosen not in b6.columns:
        raise ValueError("Column not in blast6 format.")
    b6 = b6.groupby("query")[[chosen]].max().reset_index()
    
    
    if chosen == "identity":
        b6[chosen] = round(b6[chosen]/100, 3)
    return b6


def usearch_main(usearch_result,out_name,input_file):
    b6_csv = process_blast6(usearch_result, 'identity')
    b6 = b6_csv
    identity_dict={}
    b6numpy=b6.to_numpy()
    from Bio import SeqIO
    for i in range(len(b6numpy)):
        identity_dict[b6numpy[i][0]]=b6numpy[i][1]
        
    with open('./features/biofea/'+out_name,"w",newline='') as fil1:
        fil1.write(""+","+"identity"+"\n")
        for seq_record in SeqIO.parse('./sequence/'+str(input_file),"fasta"):
            fil1.write(str(seq_record.id)+",")
            if seq_record.id not in identity_dict:
                fil1.write(str(0)+"\n")
            else:
                fil1.write(str(identity_dict[seq_record.id])+"\n")
                


# tmhmm feature
def tmhmm_fea():
    with open ('./features/biofea/tmhmmfea.csv',"w", newline="") as fil1:
        fil1.write("seq_id"+","+"AAs_TMHs"+","+"60_AAs_TMHs"+","+"N_in_cytop"+"\n")
        for line in open('./features/biofea/testing.txt'):
            if line.split(" ")[0]=="#":
                if line.split(": ")[0].split(" ")[-1]=="Length":
                    fil1.write(line.strip("\n").split(" ")[1]+",")  
                elif line.split(": ")[0].split(" ")[-2]=="in":
                    fil1.write(line.strip("\n").split(": ")[-1].strip(" ")+",")
                elif line.split(": ")[0].split(" ")[-1]=="AAs":
                    fil1.write(line.strip("\n").split(": ")[-1].strip(" ")+",")
                elif line.split(": ")[0].split(" ")[-1]=="N-in":
                    fil1.write(line.strip("\n").split(": ")[-1].strip(" ")+"\n")

# feature_combine_19
def feature_merge(feature_list,save_name):
    import numpy as np
    import pandas as pd
    n=0
    for fea in feature_list:
        #print(fea)
        n=n+1
        if n==1:           
            dfx=pd.read_csv(r'./features/biofea/'+str(fea),sep=',',header=0,index_col=0)
            dfx.pop("length")
            #print((dfx.shape[1]))
        else:
            dfn=pd.read_csv(r'./features/biofea/'+str(fea),sep=',',header=0,index_col=0)
            dfx=pd.concat([dfx,dfn],axis=1)
        dfx1=pd.DataFrame(dfx) 
        dfx1.to_csv('./features/'+str(save_name)+".csv",sep=",")
        
def biopython_features_processing(input_file):

    biopython_12(input_file)
    
    scale_main(input_file)
    
    usearch_main('testing_pdb.b6','identity_testing.csv',input_file)
    tmhmm_fea()

    feature_list=['biofea.csv', '43.csv',  '38.csv','14.csv','tmhmmfea.csv',  'identity_testing.csv',]

    feature_merge(feature_list,"19fea")

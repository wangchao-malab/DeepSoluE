# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 22:40:49 2020

@author: Administer
"""

def feature_merge_523(feature_list, save_name):
    import numpy as np
    import pandas as pd
    n=0
    for fea in feature_list:
        #print(fea)
        n=n+1
        if n==1:           
            dfx=pd.read_csv(r'./features/'+str(fea),sep=',',header=None,index_col=0)
            #print((dfx.shape[1]))
        else:
            dfn=pd.read_csv(r'./features/'+str(fea),sep=',',header=None,index_col=0)
            dfx=pd.concat([dfx,dfn],axis=1)
        dfx1=pd.DataFrame(dfx) 
        dfx1.to_csv('./features/'+str(save_name)+".csv",sep=",")

def feature_merge_119(feature_list, save_name):
    import numpy as np
    import pandas as pd
    n=0
    for fea in feature_list:
        #print(fea)
        n=n+1
        if n==1:           
            dfx=pd.read_csv(r'./features/'+str(fea),sep=',',header=0,index_col=0)
            #print((dfx.shape[1]))
        else:
            dfn=pd.read_csv(r'./features/'+str(fea),sep=',',header=0,index_col=0)
            dfx=pd.concat([dfx,dfn],axis=1)
        dfx1=pd.DataFrame(dfx) 
        dfx1.to_csv('./features/'+str(save_name)+".csv",sep=",")

def feature_merge_219(feature_list, save_name):
    import numpy as np
    import pandas as pd
    n=0
    for fea in feature_list:
        #print(fea)
        n=n+1
        if n==1:           
            dfx=pd.read_csv(r'./features/'+str(fea),sep=',',header=0,index_col=0)
            #print((dfx.shape[1]))
            index=dfx.index.to_numpy()
            #print(index)
        else:
            dfn=pd.read_csv(r'./features/'+str(fea),sep=',',header=0,index_col=0)
            dfn=dfn.set_index(index)
            #print(dfn)
            dfx=pd.concat([dfx,dfn],axis=1)
        dfx1=pd.DataFrame(dfx) 
        dfx1.to_csv('./features/'+str(save_name)+".csv",sep=",")


def phyche_merge():
    import pandas as pd
    feature_list=['AAC.csv', 'APAAC.csv', 'CTDC.csv', 'DPC.csv', 'QSOrder.csv']
    feature_merge_523(feature_list,'feature_523')

    x_523=pd.read_csv('./features/feature_523.csv',header=0,index_col=0)
    index=x_523.index
    ga_100=pd.read_csv('./model/genetic_algorithm_100_features.csv',header=0,index_col=0).to_numpy()[0]
    #print(ga_100)
    x_100=x_523.iloc[:,ga_100.tolist()].to_numpy()
    columns=["523_"+str(x) for x in ga_100]
    pd.DataFrame(x_100,index=index,columns=columns).to_csv(r'./features/feature_phyche.csv')


def feature_119():
    feature_list=['19fea.csv','feature_phyche.csv',]
    feature_merge_119(feature_list,'119feature')
    
def feature_119_minmax_scaler():
    import joblib
    import pandas as pd
    from sklearn import preprocessing
    min_max_scaler = joblib.load('./model/scaler.gz')
    dfx=pd.read_csv('./features/119feature.csv',sep=',',index_col=0,)
    dfx_ind=dfx.index
    dfx_col=dfx.columns
    scaler = min_max_scaler.transform(dfx)
    scaler=pd.DataFrame(scaler,index=dfx_ind,columns=dfx_col)
    scaler.to_csv(r'./features/119feature_scaler.csv',sep=',') 

def feature_219():
    feature_list=['119feature_scaler.csv','w2v_feature.csv',]
    feature_merge_219(feature_list,'219feature')
        
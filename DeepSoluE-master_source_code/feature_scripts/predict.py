# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 22:40:49 2020

@author: Administer
"""
def result_prediction(output_file):
    import numpy as np
    import pandas as pd
    import os
    import sys
    import tensorflow as tf
    from tensorflow import keras
    os.environ["TF_CPP_MIN_LOG_LEVEL"]='3'
    feature_length=219
    results_id_proba_site=[]
    
    pred_proba=[]
    x_test=pd.read_csv('./features/219feature.csv',sep=',',header=0,index_col=0)
    seq_id=list(x_test.index)
    x_test=x_test.to_numpy()
    results_id_proba_site.append(seq_id)
    x_test_length=len(x_test)
    x_test=np.float32(x_test.reshape(x_test_length, feature_length,1))
    for model_num in range(1,11,1):
        tf_model = keras.models.load_model("./model/tf_model/model_"+str(model_num)+".hdf5")        
        predict_proba=tf_model.predict(x_test)[:,1]
        pred_proba.append(list(predict_proba))
    
    pred_proba_ave=np.mean(np.array(pred_proba).T,axis=1)
    pred_proba_ave=[np.round(x,4) for x in pred_proba_ave]
    results_id_proba_site.append(pred_proba_ave)    
      
    pred_class=np.array(pred_proba).T
    sss=pred_class.reshape(x_test_length*10)
    sss=list(map(lambda x: 1 if x>0.4 else 0,sss)) 
    sss=np.array(sss).reshape(x_test_length,10)
    pred_class=np.sum(np.array(sss),axis=1)    
    pred_class=list(map(lambda x: "soluble" if x>5 else "insoluble",pred_class))
    results_id_proba_site.append(pred_class)
    pd.DataFrame(np.array(results_id_proba_site).T,columns=["sequence_id","predicted_probability","result"]).to_csv('./results/'+output_file)
    
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 22:40:49 2020

@author: Administer
"""

def DeepSoluE(input_file,output_file):
    import os,sys
    work_path=os.path.abspath('.')
    os.chdir(work_path)
    sys.path.append('./feature_scripts/')
    # print(os.getcwd())
    import AAC, APAAC, CTDC, DPC, QSOrder, Biofea, w2v_kmer_corpus_feature
    import sequence_read_save
    import feature_combine, predict
    import tmhmm_usearch
    print("step_1: sequence checking......")
    input_file=input_file
    fastas = sequence_read_save.read_protein_sequences(input_file)

    
    print("step_2: sequence encoding......")
    sequence_read_save.file_remove()
    QSOrder.QSOrder_feature(fastas)
    AAC.AAC_feature(fastas)
    APAAC.APAAC_feature(fastas)
    CTDC.CTDC_feature(fastas)
    DPC.DPC_feature(fastas)
    tmhmm_usearch.usearch(input_file)
    tmhmm_usearch.tmhmm(input_file)
    Biofea.biopython_features_processing(input_file)
    feature_combine.phyche_merge()
    w2v_kmer_corpus_feature.w2v_kmer_corpus(input_file)
    w2v_kmer_corpus_feature.w2v_features()
    feature_combine.feature_119()
    feature_combine.feature_119_minmax_scaler()
    feature_combine.feature_219()
    print("step_3: result preiction.......")
    predict.result_prediction(output_file)
    sequence_read_save.file_remove()
if __name__=='__main__':
    import argparse
    import pandas as pd
    import os
    work_path=os.path.abspath('.')
    os.chdir(work_path)
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inputfile", help="input fasta file",)
    parser.add_argument("-o","--outputfile", help="output fasta file",)
    args = parser.parse_args()
    input_file=args.inputfile
    output_file=args.outputfile
    print("work launched")
    DeepSoluE(input_file,output_file)
    print("job finished !")
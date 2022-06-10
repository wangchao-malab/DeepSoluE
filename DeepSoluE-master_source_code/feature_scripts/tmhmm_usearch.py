# -*- coding: utf-8 -*-

def tmhmm(input_file):
	import os
	#print("working dir: ", os.getcwd()) 
	#os.chdir('/home/wangchao/web_server/nodejs/DeepSoluE-master/softs/tmhmm-2.0c')
	comm_line="./softs/tmhmm-2.0c/bin/tmhmm  ./sequence/"+str(input_file)+" >>./features/biofea/testing.txt"
	#print(comm_line)
	os.system(comm_line)


def usearch(input_file):
	import os
	#print("working dir: ", os.getcwd()) 
	comm_line="./softs/usearch/usearch -usearch_global ./sequence/"+str(input_file)+"  -db ./softs/usearch/Ecoli_xray_nmr_pdb_no_nesg_simple_id.fasta  -id 0 -blast6out ./features/biofea/testing_pdb.b6"
	#print(comm_line)
	os.system(comm_line)
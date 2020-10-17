#script for evaluate the performance of cross_validation procedures 
import performance as p
import os, re
import sys, math, statistics
import numpy as np
if __name__=='__main__':
	filename1=sys.argv[1]
	filename2=sys.argv[2]
	filename3=sys.argv[3]
	filename4=sys.argv[4]
	filename5=sys.argv[5]
	mainresults1=p.conf_mat(filename1)
#	print(mainresults1)
	mainresults2=p.conf_mat(filename2)
	mainresults3=p.conf_mat(filename3)
	mainresults4=p.conf_mat(filename4)
	mainresults5=p.conf_mat(filename5)
#	(SOV_H_average, SOV_E_average, SOV_C_average, q3, sen_H, ppv_H, sen_E, ppv_E, sen_C, ppv_C, mcc_H, mcc_E, mcc_C )
	listt=['SOV_H_average', 'SOV_E_average', 'SOV_C_average', 'q3', 'sen_H', 'ppv_H', 'sen_E', 'ppv_E', 'sen_C', 'ppv_C', 'mcc_H', 'mcc_E', 'mcc_C']
	for i in range(0,13):
		average= (mainresults1[i]+mainresults2[i]+mainresults3[i]+ mainresults4[i]+ mainresults5[i])/5
		sd=statistics.stdev([mainresults1[i],mainresults2[i],mainresults3[i], mainresults4[i], mainresults5[i]])
		print('mean_cv %s:'%(listt[i]),average)
		print('SE_cv %s:' %(listt[i]), sd/(math.sqrt(5)))

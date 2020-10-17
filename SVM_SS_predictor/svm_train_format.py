##svm script for creating the input format for svm-train program
import os 
import numpy as np
from sklearn import svm

def make_XY(train_id):
	d={'H':1, 'E':2, '-':3}
	with open('/datacomplete.txt' , 'r') as file2:
		file2_list= file2.readlines()
		file2_list_strip=[i.strip() for i in file2_list]
		with open(train_id, 'r') as idtext:
			for line in idtext:
				idname= line.rstrip()
				if idname + '.norm.prof' in os.listdir('/seq_profile_norm/'):
					with open('/seq_profile_norm/'+ idname + '.norm.prof') as matrix:
						matrixline= matrix.readlines()
						L=[]
						for i in range(1, len(matrixline)):
							line=matrixline[i].split()[2:]
							num=[float(ch) for ch in line]
							L.append(num)
						A=np.array(L)
						if np.sum(A)!=0:
							ss_string=file2_list_strip[file2_list_strip.index('>'+ idname)+ 1]
							counter=-1
							An=np.vstack((np.zeros((8,20)),A,np.zeros((8,20))))
							for i in range(0,A.shape[0]):
								W=An[i:i+17]
								counter+=1
								if np.sum(W)!=0:
									W_string=''
									W_merged=''
									W_string+='{0} '.format(d[ss_string[counter]])
									for i in W:
										for j in i:
											W_merged+="{0} ".format(j)
									for i in range(0,340):	
										if W_merged.split(" ")[i]!='0.0':
											W_string+="{0}:{1} ".format(i, W_merged.split(" ")[i])
									print(W_string)

	
if __name__=="__main__":
	myfile= os.sys.argv[1]
	make_XY(myfile)

                                      
								

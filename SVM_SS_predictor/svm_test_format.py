##SVM script for creating the input format for svm-predict
#!/usr/bin/env python2
import os 
import numpy as np

def make_XY(pred_id):
	d={'H':1, 'E':2, '-':3}
	with open('blind_data_complete.fasta' , 'r') as file2:
		file2_list= file2.readlines()
		file2_list_strip=[i.strip() for i in file2_list]
		with open(pred_id, 'r') as idtext:
			for line in idtext:
				idname= line.rstrip().split('.')[0]
				if idname + '.norm.prof' in os.listdir('/seq_profile_norm/'):
					with open('/seq_profile_norm/'+ idname + '.norm.prof') as matrix:
						matrixline= matrix.readlines()
						L=[]
						for i in range(1, len(matrixline)):
							line=matrixline[i].split()[2:]
							num=[float(ch) for ch in line]
							L.append(num)
						A=np.array(L)
                                                #one-hot encode technique
						if np.sum(A)==0:
							alphabet='ARNDCQEGHILKMFPSTWYVX'
							char_to_int=dict((c,i) for i,c in enumerate(alphabet))
							fastaseq=file2_list_strip[file2_list_strip.index('>'+ idname)+ 1]
							data_seq=fastaseq
							integer_encoded=[char_to_int[char] for char in data_seq]
							onehot_encoded = list()
							for value in integer_encoded:
								zero_matr = [0 for _ in range(len(alphabet)-1)]
								if value !=20:
									zero_matr[value] = 1
								onehot_encoded.append(zero_matr)
								A=np.array(onehot_encoded)

						ss_string=file2_list_strip[file2_list_strip.index('>'+ idname)+ 2]

						counter=-1
						An=np.vstack((np.zeros((8,20)),A,np.zeros((8,20))))
						final_string=''
						for i in range(0,A.shape[0]):
							W=An[i:i+17]
							counter+=1

							W_string=''		
							W_merged=''
							W_string+='{0} '.format(d[ss_string[counter]])
							for i in W:
								for j in i:
									W_merged+="{0} ".format(j)
							for i in range(0,340):
								W_string+="{0}:{1} ".format(i, W_merged.split(" ")[i])
							final_string= final_string + W_string + '\n'
						with open(idname + '.SVMtest_blind.dat', 'w') as pp:
							pp.write(final_string)

	
if __name__=="__main__":
	myid= os.sys.argv[1]
	make_XY(myid)

                                      
								

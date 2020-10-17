##GOR predict script for predicting blind set
import os
import numpy as np
import gor_train as gr

def gor_pred(profdir,fastadir):
	with open('gor_blind_result.txt', 'w') as filefin:
		with open('/home/um79/project/blind_test/blind_data_complete.fasta' , 'r') as file2:
			file2_list= file2.readlines()
			file2_list_strip=[i.strip() for i in file2_list]
			for filepssm in os.listdir(profdir):
				idname=os.path.splitext(filepssm)[0].split('.')[0]
				seq=open(os.path.join(fastadir,idname + '.fa'))
				fastaseq= seq.readlines()[1].rstrip()[1:]
				ss_seq=file2_list_strip[file2_list_strip.index('>'+ idname)+ 2]
				with open(os.path.join(profdir,filepssm)) as matrix:
					matrixline=matrix.readlines()
					L=[]
					for i in range(1, len(matrixline)):
						line=matrixline[i].split()[2:]
						num=[float(ch) for ch in line]
						L.append(num)
					A=np.array(L)
					if np.sum(A)==0:
						alphabet='ARNDCQEGHILKMFPSTWYVX'
						char_to_int=dict((c,i) for i,c in enumerate(alphabet))
						data_seq=fastaseq
						integer_encoded=[char_to_int[char] for char in data_seq]
						onehot_encoded = list()
						for value in integer_encoded:
							zero_matr = [0 for _ in range(len(alphabet)-1)]
							if value !=20:
								zero_matr[value] = 1
							onehot_encoded.append(zero_matr)
							A=np.array(onehot_encoded)
					An=np.vstack((np.zeros((8,20)),A,np.zeros((8,20))))
					main_matrices=gr.gor('/home/um79/project/train_fasta/seq_profile_norm/')
					H=main_matrices[0]
					E=main_matrices[1]
					C=main_matrices[2]
					R=main_matrices[3]
					stru=main_matrices[4] #list with helix,coil,strand
					ss_string=''
					for i in range(0, A.shape[0]):
						W= An[i:i+17]
						H_score= np.true_divide(np.multiply(W,H),R)
						H_score2=np.true_divide(H_score, stru[0])
						H_sum=np.sum(H_score2)
						C_score= np.true_divide(np.multiply(W,C),R)
						C_score2=np.true_divide(C_score, stru[1])
						C_sum=np.sum(C_score2)
						E_score= np.true_divide(np.multiply(W,E),R)
						E_score2=np.true_divide(E_score, stru[2])
						E_sum=np.sum(E_score2)				
						stru_list=[H_sum, C_sum, E_sum]
						max_index= stru_list.index(max(stru_list))
						if max_index ==0:
							ss_string+= 'H'
						if max_index==1:
							ss_string+= '-'
						if max_index==2:
							ss_string+= 'E'
					filefin.write('\n'.join(['>'+idname, ss_string, ss_seq]) + '\n')

if __name__== '__main__':
	profdirect= os.sys.argv[1]
	fastadirect= os.sys.argv[2]
	gor_pred(profdirect,fastadirect)

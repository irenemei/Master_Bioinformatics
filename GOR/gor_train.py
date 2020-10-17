##GOR train script for training the Jpred4 set
import os
import numpy as np

def gor(profdir):
	count=0
	with open('datacomplete.txt' , 'r') as file2:
		file2_list= file2.readlines()
		file2_list_strip=[i.strip() for i in file2_list]
		H=np.zeros((17,20))
		E=np.zeros((17,20))
		C=np.zeros((17,20))
		R=np.zeros((17,20))
		listruct=[0,0,0]
		for  file1 in os.listdir(profdir):
			idname=os.path.splitext(file1)[0].split('.')[0]
			ss_string=file2_list_strip[file2_list_strip.index('>'+ idname)+ 1]	#+1 per train datacomplete +2 per blind
			with open(os.path.join(profdir, file1)) as matrix:
				matrixline= matrix.readlines()
				L=[]
				for i in range(1, len(matrixline)):
					line=matrixline[i].split()[2:]
					num=[float(ch) for ch in line]	
					L.append(num)			
				A=np.array(L)
			
				if np.sum(A)!=0:
					An=np.vstack((np.zeros((8,20)),A,np.zeros((8,20))))
					for i in range(0,A.shape[0]):
						W=An[i:i+17]
						R=R+W
						if ss_string[i]=='H':
							H=H+W
							listruct[0]=listruct[0] + 1
						if ss_string[i]=='-':
							C=C+W
							listruct[1]=listruct[1] + 1
						if ss_string[i]=='E':
							E=E+W
							listruct[2]=listruct[2] + 1

					H_f= np.true_divide(H,A.shape[0])
					C_f= np.true_divide(C,A.shape[0])
					E_f= np.true_divide(E,A.shape[0])
					R_f= np.true_divide(R,A.shape[0])

						
		for i in H:
			print (i)
#	return(H_f,E_f,C_f,R_f,listruct)	
			
if __name__=="__main__":
        p= os.sys.argv[1]
        gor(p)

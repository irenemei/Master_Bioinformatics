## Plots statistical analysis
import sys
import os, sys
import numpy as np
import seaborn as sns; sns.set()
def make_ss_dict(datafile):
	t={'H':0, '-':0, 'E':0}
	with open(datafile, 'r') as f:
		fline= f.readlines()
		j=1
		while (j<4044):
			for i in fline[j].rstrip():				
				t[i]=t[i]+2
			j=j+3
	return t
def make_ss_pie(g):
	from pylab import figure, title, pie, savefig
	import matplotlib.pyplot as plt	
	ss=["Helix", "Coil", "Strand"]
	count= g.values()
	colors="blue", "green", "magenta"
	def get_percent(value): return "%4.1f%%" % (value)
	figure(1,figsize=(10,10))
	plt.rcParams['font.size'] = 15
	ppie=pie(count, labels=ss, colors=colors, autopct=get_percent)
	savefig('pie_chart.png', dpi=150)
	return

def make_res_ss_dict(datafile):
	
	d= {
	"A":{"H":0,"E":0,"-":0,"O":0},"R":{"H":0,"E":0,"-":0,"O":0},"N":{"H":0,"E":0,"-":0,"O":0},"D":{"H":0,"E":0,"-":0,"O":0},
	"C":{"H":0,"E":0,"-":0,"O":0},"Q":{"H":0,"E":0,"-":0,"O":0},"E":{"H":0,"E":0,"-":0,"O":0},"G":{"H":0,"E":0,"-":0,"O":0},
	"H":{"H":0,"E":0,"-":0,"O":0},"I":{"H":0,"E":0,"-":0,"O":0},"L":{"H":0,"E":0,"-":0,"O":0},"K":{"H":0,"E":0,"-":0,"O":0},
	"M":{"H":0,"E":0,"-":0,"O":0},"F":{"H":0,"E":0,"-":0,"O":0},"P":{"H":0,"E":0,"-":0,"O":0},"S":{"H":0,"E":0,"-":0,"O":0},
	"T":{"H":0,"E":0,"-":0,"O":0},"W":{"H":0,"E":0,"-":0,"O":0},"Y":{"H":0,"E":0,"-":0,"O":0},"V":{"H":0,"E":0,"-":0,"O":0}
	}	 
	with open(datafile, 'r') as f:
		fline= f.readlines()						
		j=1
		while (j<4044):
			u=0
			for char in fline[j].rstrip():
				if char!="X":
					d[char]["O"]+=1
					s=fline[j-1].rstrip()[u]
					d[char][s]+=1
					u+=1
			j=j+3
	return d

def make_histo(m):
	import numpy
	from pylab import figure, title, xlabel, ylabel, xticks, yticks, bar, legend, axis, savefig
	H = [ float(x["H"]) for x in m.values()]
	E = [ float(x["E"]) for x in m.values()]
	C = [ float(x["-"]) for x in m.values()]
	O = [ float(x["O"]) for x in m.values()]
	tot_H=sum(H)
	tot_E=sum(E)
	tot_C=sum(C)
	tot_AA=sum(O) 
	H_freq= [(h/tot_H)*100 for h in H]
	E_freq= [(e/tot_E)*100 for e in E]
	C_freq= [(c/tot_C)*100 for c in C]
	O_freq= [(o/tot_AA)*100 for o in O]
	residues = m.keys()
	figure()
	title('Residue composition')
	xlabel('Residues')
	ylabel('Residue frequency (%)')
	x1 = numpy.array(range(2,41,2))
	x2 = [x - 0.25 for x in x1]
	x3 = [x - 0.25 for x in x2]
	x4 = [x - 0.25 for x in x3]
	xticks(x2, residues)
	yticks(range(2,17,2), ["2%", "4%", "6%", "8%", "10%", "12%", "14%", "16%" ])
	bar(x1, C_freq, width=0.25, color="blue", label="Coil")
	bar(x2, E_freq, width=0.25, color="yellow", label="Strand")
	bar(x3, H_freq, width=0.25, color="red", label="Helix")
	bar(x4, O_freq, width=0.25, color="green", label="Overall")
	legend()
	axis([1.0, 43, 0, 17])
	savefig('histogram.png')
	return			
	
	
def taxonomy(tax_file):
	from pylab import figure, title, pie, savefig
	import matplotlib.pyplot as plt
	with open(tax_file, 'r') as fi:
		a=fi.readlines()
		a_strip= [i.rstrip() for i in a]
		#print a_strip
		d = {x:a_strip.count(x) for x in a_strip} #string of the dictionary is without \n			
		d_filter=dict()
		others=0
		# Iterate over all the items in dictionary and filter items which has even keys
		for (key, value) in d.items():
		   # Check if key is even then add pair to new dictionary
			if value >= 40:
				if value != 87 :
					d_filter[key] = value
				else:
					others= others + value
			else:
				others= others + value
		d_filter["Other"]= others
		#print d_filter		
	def get_percent(value): return "%4.1f%%" % (value)
	figure(2, figsize=(15,15))
	colors=["blue","red", "yellow", "orange", "violet", "green", "magenta"]
	#plt.axis["Equal"]
	print (d_filter)
	a=pie(d_filter.values(), labels= d_filter.keys(), colors=colors, autopct=get_percent)
	title("Taxonomic classification: species ")
	savefig('pie_taxonomy.png', dpi=150)
	
def kingdom_pie():
	from pylab import figure, title, pie, savefig
	import matplotlib.pyplot as plt	
	taxon=["Bacteria","Eukaryota", "Archaea", "Viruses", "Other"]
	count=[60,70,7,11,2]
	def get_percent(value):
		        return "%4.1f%%" % (value)
	figure((3), figsize=(15,15))
	plt.rcParams['font.size'] = 21
	pie(count, labels=taxon, autopct=get_percent)
	title('Taxonomic classification: kingdom')	
	savefig('pie_kingdom.png', dpi=150)


def scop_pie():
	from pylab import figure, title, pie, savefig
	import matplotlib.pyplot as plt
	structure= ["Alpha + beta (a+b)"," All alpha","All beta","Alpha and beta (a/b)","Small proteins", "Multi-domain","Other"]
	count=[323,302,218,205,42,31,24]
	def get_percent(value):
		        return "%4.1f%%" % (value)
	figure((4),figsize=(15,15))
	plt.rcParams['font.size'] = 21
	pie(count, labels=structure,autopct=get_percent )
	title('SCOP classification')	
	savefig('pie_SCOP.png', dpi=150)
def heatmap(): 
	with open('/home/um79/project/Statics/datacomplete.txt' , 'r') as file2:
		file2_list= file2.readlines()
		file2_list_strip=[d.strip() for d in file2_list]
		H=np.zeros((21,17))
		E=np.zeros((21,17))
		C=np.zeros((21,17))
		j=1
		alphabet='ARNDCQEGHILKMFPSTWYVX'
		char_to_int=dict((c,a) for a,c in enumerate(alphabet))          #con X in 20
		while (j<4044):
			fa_seq=file2_list_strip[j+1]
			ss_seq=file2_list_strip[j]
			integer_encoded =[char_to_int[char]for char in fa_seq] 
			W_index=0
			while(W_index+17< len(fa_seq)):
				W=integer_encoded[W_index:W_index+17]
				W_ss=ss_seq[W_index:W_index+17]
				for i in range(0,17):
					if W_ss[i]=='H':
						H[W[i],i]+=1
					elif W_ss[i]=='E':
						E[W[i],i]+=1
					elif W_ss[i]=='-':
						C[W[i],i]+=1
				W_index+=1
			j=j+3
		h_heatmap = sns.heatmap(H)
		print(h_heatmap)
		e_heatmap = sns.heatmap(E)
		c_heatmap = sns.heatmap(C)

 
			
if __name__=="__main__":
	data_file= sys.argv[1]
	source_file= sys.argv[2] 	#file with taxonomy class related to organism of the protein id	
	make_ss_pie(a)
	dict_res=make_res_ss_dict(data_file)
	make_histo(dict_res)
	taxonomy(source_file) 
	kingdom_pie()
	scop_pie()
	heatmap()


import os, re
import sys, math
import numpy as np

def conf_mat(filename):
	cm=np.array([[0.0,0.0,0.0],[0.0,0.0,0.0], [0.0,0.0,0.0]])
	with open (filename) as fn:
		fn_lines= fn.readlines() 	#lista di righe
		fn_lines=[i.strip() for i in fn_lines]
		i=1
		count=0
		sovH_numerator=0
		counter_0_sov_H=0
		sovE_numerator=0
		counter_0_sov_E=0
		sovC_numerator=0
		counter_0_sov_C=0
		while i<len(fn_lines):
			o_ss=fn_lines[i]
			p_ss=fn_lines[i+1]
			N=len(o_ss)
			for t in range(0, len(o_ss)):
				if o_ss[t]==p_ss[t]:
					if o_ss[t]=='H': cm[0][0]+=1
					elif o_ss[t]=='E': cm[1][1]+=1
					elif o_ss[t]=='-': cm[2][2]+=1
				elif o_ss[t]=='H' and p_ss[t]=='E': cm[1][0]+=1
				elif o_ss[t]=='H' and p_ss[t]=='-': cm[2][0]+=1
				elif o_ss[t]=='E' and p_ss[t]=='H': cm[0][1]+=1
				elif o_ss[t]=='E' and p_ss[t]=='-': cm[2][1]+=1
				elif o_ss[t]=='-' and p_ss[t]=='H': cm[0][2]+=1
				elif o_ss[t]=='-' and p_ss[t]=='E': cm[1][2]+=1
			p=re.compile('H+')
			So = p.finditer(o_ss)	
			Sp = p.finditer(p_ss)
			o_list_span=[h.span() for h in So]
			o_list_range=[range(a,b) for (a,b) in o_list_span]
			p_list_span=[h.span() for h in Sp]
			p_list_range=[range(a,b) for (a,b) in p_list_span]
			counter=0	#count minov*2l/maxov since in formula is a summation
			counter_H=0
			for h in o_list_range:
				a=set(h)
				counter_H += len(a) #is a counter for the H
				for j in p_list_range:
					b=set(j)
					if list(a&b)!=[]:
						c=list(a&b)
						minov=len(c)	#numero effettivo di elementi che overlappano

						maxov=(len(a)+len(b)-minov) #ex. (2,8)
						l=len(a)
						l2=len(b)
						delta=min((maxov-minov), minov,l/2,l2/2)
						counter+=((minov+delta)*l)/maxov
			if counter_H!=0:
				SOV_H=100*(1/counter_H)*counter
				sovH_numerator+= SOV_H
			else:
				counter_0_sov_H+= 1				
#####SOV_E#####
			p=re.compile('E+')
			SoE = p.finditer(o_ss)   #find all occurences in the string but save also other attributes that can be usefull after like match.span, match.group
			SpE = p.finditer(p_ss)

			o_list_span_E=[e.span() for e in SoE]
			o_list_range_E=[range(a,b) for (a,b) in o_list_span_E]
			#print('predicted')
			p_list_span_E=[e.span() for e in SpE]
			p_list_range_E=[range(a,b) for (a,b) in p_list_span_E]
			counter2=0
			counter_E=0
			for e in o_list_range_E:
				a2=set(e)
				counter_E += len(a2) #is a counter for the H
				for j in p_list_range_E:
					b2=set(j)
					if list(a2&b2)!=[]:
						c2=list(a2&b2)
						minov=len(c2)    #numero effettivo di elementi che overlappano
						maxov=(len(a2)+len(b2)-minov) #ex. (2,8)
						l=len(a2)
						l2=len(b2)
						delta=min((maxov-minov), minov,l/2,l2/2)
						counter2+=((minov+delta)*l)/maxov
			if counter_E!=0:	
                                SOV_E=100*(1/counter_E)*counter2
                                sovE_numerator+= SOV_E
			else:
                                counter_0_sov_E+= 1
			#	print(match2.span(), match2.group())		
##SOVC##
			p=re.compile('-+')
			SoC = p.finditer(o_ss)   #find all occurences in the string but save also other attributes that can be usefull after like match.span, match.group
			SpC = p.finditer(p_ss)
			#print (match.span(),match.group())
			o_list_spanC=[c.span() for c in SoC] 
			o_list_rangeC=[range(a,b) for (a,b) in o_list_spanC]
			p_list_spanC=[c.span() for c in SpC]
			p_list_rangeC=[range(a,b) for (a,b) in p_list_spanC]
			counter3=0
			counter_C=0
			for c in o_list_rangeC:
				a3=set(c)
				counter_C += len(a3)
				for j in p_list_rangeC:
					b3=set(j)
					if list(a3&b3)!=[]:
						c3=list(a3&b3)
						minov=len(c3)    

						maxov=(len(a3)+len(b3)-minov) 
						l=len(a3)
						l2=len(b3)
						delta=min((maxov-minov), minov,l/2,l2/2)
						counter3+=((minov+delta)*l)/maxov
			if counter_C!=0:
                                SOV_C=100*(1/counter_C)*counter3
                                sovC_numerator+= SOV_C
			else:
                                counter_0_sov_C+= 1
			i+=3
		SOV_H_average= sovH_numerator /((len(fn_lines)/3)-counter_0_sov_H)
		SOV_E_average= sovE_numerator /((len(fn_lines)/3)-counter_0_sov_E)
		SOV_C_average= sovC_numerator /((len(fn_lines)/3)-counter_0_sov_C)
	
			 
		

	q3=np.trace(cm)/np.sum(cm)
	sen_H= cm[0][0]/np.sum(cm[:,0])	
	ppv_H= cm[0][0]/np.sum(cm[0,:])
	sen_E= cm[1][1]/np.sum(cm[:,1])
	ppv_E= cm[1][1]/np.sum(cm[1,:])
	sen_C= cm[2][2]/np.sum(cm[:,2])
	ppv_C= cm[2][2]/np.sum(cm[2,:])

	cH=cm[0,0]
	nH=cm[1:3,1:3].sum()
	uH=cm[1:3,0].sum()
	oH=cm[0,1:3].sum()
	mcc_H= ((cH*nH)-(oH*uH))/math.sqrt((cH+oH)*(cH+uH)*(nH+oH)*(nH+uH))
	cE=cm[1,1]
	nE=(cm[:,0].sum() - cm[1,0])+ (cm[:,2].sum() - cm[1,2])
	uE=cm[:,1].sum() - cm[1,1]
	oE=cm[1,:].sum() - cm[1,1]
	mcc_E= ((cE*nE)-(oE*uE))/math.sqrt((cE+oE)*(cE+uE)*(nE+oE)*(nE+uE))
	cC=cm[2,2]
	nC=cm[0:2,0:2].sum()
	uC=cm[0:2,2].sum()
	oC=cm[2,0:2].sum()
	mcc_C= ((cC*nC)-(oC*uC))/math.sqrt((cC+oC)*(cC+uC)*(nC+oC)*(nC+uC))
#			print (count)
	print('SOV_H:', SOV_H_average, 'SOV_E:', SOV_E_average,'SOV_C:', SOV_C_average,'acc:', q3, 'sen_H:',sen_H,'ppv_H:', ppv_H, 'sen_E:', sen_E, 'ppv_E:',ppv_E, 'sen_C:', sen_C, 'ppv_C:', ppv_C, 'mcc_H:', mcc_H,'mcc_E:', mcc_E,'mcc_C:', mcc_C )
#	return(SOV_H_average, SOV_E_average, SOV_C_average, q3, sen_H, ppv_H, sen_E, ppv_E, sen_C, ppv_C, mcc_H, mcc_E, mcc_C )	##return is used for the performance_cv.py	
if __name__=='__main__':
	inp=sys.argv[1]
	conf_mat(inp)

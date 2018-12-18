import numpy as np
import pandas as pd
from scipy.constants import codata
from RF_mathematica import *

hc=codata.value('inverse meter-electron volt relationship')*1E+2 #eV.cm

def iter_calc(l2, react_S_B, react_L_B, I_B, term, J, L, S, level, skipped_level, second_ce, df_double_ce, ele_B, ele_A, prod_S_A, react_S_A, react_L_A, charge_state, T_B, dist):

	ele_sym_B=ele_B.symbol
	ele_sym_A=ele_A.symbol
	I_A=ele_A.ionenergies[1]
	m_B=ele_B.mass

	try:
		if second_ce:
			raise ValueError("second_ce true")

		levels_pops_detunes_f = open('results/Z'+str(ele_B.atomic_number)+'/levels_pops_detunes'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist),'rb')
		levels_pops_detunes=np.array(list(np.loadtxt(levels_pops_detunes_f, delimiter=';',dtype=float)))

		levels_cs_f = open('results/Z'+str(ele_B.atomic_number)+'/levels_cs'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist),'rb')
		levels_cs=np.array(list(np.loadtxt(levels_cs_f, delimiter=';',dtype=float)))

		levels_pops_detunes_initial=levels_pops_detunes

		return 0.0, levels_pops_detunes, levels_cs, df_double_ce

	except Exception as exception:
		print(exception)
		#all of  below
		###ction for each level####
		cross_sections=[]
		energy_levels=[]
		bare_cross_sections=[]

		# print("for index,x in enumerate(level[0:]):")
		# print(level[0:])

		for index,x in enumerate(level[0:]):
			# print(x)
			x=x.replace("[","").replace("]","").replace("?","").replace(" ","")
			I_B_ex=I_B-float(x)*hc #eV

			try:
				cs=cross_section(I_A,I_B_ex,m_B,T_B)
			except Exception as e:
				print(e)
				cs=0.0
				print("cs=0.0 going to next level")

			prod_S_B=S[index]
			prod_L_B=L[index]

			if (prod_L_B or prod_S_B) == '*':
				print("*",prod_L_B,prod_S_B)
				continue

			#products/reactant multiplicity
			f_s=((2*prod_S_B+1)*(2*prod_S_A+1))/((2*react_S_B+1)*(2*react_S_A+1))
			if f_s > 1:
				f_s=1.0
			f_l=((2*prod_L_B+1)*(2*prod_S_A+1))/((2*react_L_B+1)*(2*react_L_A+1))
			if f_l > 1:
				f_l=1.0

			f=f_s*f_l

			print(f,x,cs)
			energy_levels.append(float(x)) # add 0.001 if lower level energy already there?
			if cs == None:
				cs=0.0
			cross_sections.append(f*cs)
			bare_cross_sections.append(cs) # saves bare ctions for the resonance plot later..

		###cross-section for each level####

		# print(cross_sections)
		print("total sqrt cross-section:", sum(cross_sections), "(cm)")

		if not second_ce:
			np.savetxt('results/Z'+str(ele_B.atomic_number)+'/total_cross_section'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist), np.array([sum(cross_sections)]), delimiter=';')

			cs_norm = [float(i)/sum(cross_sections) for i in cross_sections]

			# levels_pops_detunes=[] # [level,[[pop1,detune1],[pop1,detune1]]] structure or [[level1,pop1,detune1],[level1,pop1,detune1]]]
			levels_pops_detunes=np.array(list(zip(energy_levels,cs_norm,[0.0] * len(energy_levels))), dtype=float)
			levels_cs=np.array(list(zip(energy_levels, bare_cross_sections)), dtype=float)
			levels_pops_detunes_initial=levels_pops_detunes
			print("initial states populated")
			print(levels_pops_detunes)

			print("total cross-section:", sum(cross_sections), "(cm^2)")

			np.savetxt('results/Z'+str(ele_B.atomic_number)+'/levels_pops_detunes'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist), levels_pops_detunes_initial, delimiter=';')
			np.savetxt('results/Z'+str(ele_B.atomic_number)+'/levels_cs'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist), levels_cs, delimiter=';')

		else:
			cs_norm = [float(i)/sum(cross_sections) for i in cross_sections]
			levels_pops_detunes=np.array(list(zip(energy_levels,cs_norm,[0.0] * len(energy_levels))), dtype=float)
			levels_cs=np.array(list(zip(energy_levels, bare_cross_sections)), dtype=float)
			levels_pops_detunes_initial=levels_pops_detunes

			d_dict={"IIlevel":[l2]*len(levels_pops_detunes[:,0]),"Ilevel":levels_pops_detunes[:,0], "normpop":levels_pops_detunes[:,1], "cs":levels_cs[:,1]}
			# d_dict={"level":levels_pops_detunes[:,0], "pop":levels_pops_detunes[:,1]}
			df_this_I_B = pd.DataFrame(data=d_dict)
			df_double_ce=pd.concat([df_double_ce, df_this_I_B])

			file_string='results/Z'+str(ele_B.atomic_number)+'/second_ce/levels_pops_detunes'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist)
			df_double_ce.to_csv(file_string)
			print(df_double_ce)

		return sum(cross_sections), levels_pops_detunes, levels_cs, df_double_ce

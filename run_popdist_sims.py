import numpy as np
import pandas as pd
from mendeleev import element
from scipy.constants import codata

from read_database import load_levels_lines
from RF_mathematica import cross_section
from cumsum_diff import cumsum_diff
from iter_calc import *

import itertools
import matplotlib.pyplot as plt, mpld3

import scipy as sp
import sympy as sympy
from sys import *
import os
import ast
import scipy.stats
import matplotlib.cm as cm

sim_list='sim_list_all.csv'

atomionspins_file = open("datatables/nistatomionspins.csv",'rb')
atomionspins = np.genfromtxt(atomionspins_file,delimiter = '\t',dtype=str,skip_header=1,autostrip=1)


def popsim(ele_sym_A,ele_sym_B,T_B,dist,sidepeaks_transition,sidepeaks_collinear, time_steps, database, charge_state, double_ce):

	print(sidepeaks_transition)
	sidepeaks_sim=bool(sidepeaks_transition)
	print(sidepeaks_sim)

	try:
		ele_num_B=int(ele_sym_B)
		ele_sym_B=element(ele_num_B).symbol

	except:
		pass

	ele_A=element(ele_sym_A) #e.g. Na, K, Li
	ele_B=element(ele_sym_B)

	if not os.path.exists("results/Z"+str(ele_B.atomic_number)):
		os.makedirs("results/Z"+str(ele_B.atomic_number))

	row = np.where(atomionspins[:,1] == ele_sym_A)

	react_L_A=float(atomionspins[:,2][row][0]) #I
	react_S_A=float(atomionspins[:,3][row][0])
	prod_L_A=float(atomionspins[:,4][row][0]) # II
	prod_S_A=float(atomionspins[:,5][row][0])

	row = np.where(atomionspins[:,1] == ele_sym_B)

	react_L_B=float(atomionspins[:,4][row][0]) # II
	react_S_B=float(atomionspins[:,5][row][0])

	show_notrans=0

	m_B=ele_B.mass
	I_A=ele_A.ionenergies[1]

	B_string=ele_sym_B+'I'*charge_state

	if dist !=0:
		I_B, term, J, L, S, level, skipped_level,  E_k_probs = load_levels_lines(B_string, 0, dist, database, ele_B, charge_state, I_A)
	else:
		I_B, term, J, L, S, level, skipped_level = load_levels_lines(B_string, 0, dist, database, ele_B, charge_state, I_A)

	print(B_string+ " levels and lines loaded")


	df_double_ce=pd.DataFrame()

	tot_cs, levels_pops_detunes, levels_cs, df_double_ce = iter_calc([0.0], react_S_B, react_L_B, I_B, term, J, L, S, level, skipped_level, 0, df_double_ce, ele_B, ele_A, prod_S_A, react_S_A, react_L_A, charge_state, T_B, dist)

	levels_pops_detunes_initial = levels_pops_detunes

	I_B2, term2, J2, L2, S2, level2, skipped_level2 = I_B, term, J, L, S, level, skipped_level #save charge state 2 for iteration

	second_ce_levels, second_ce_cs=[], []
	if double_ce:
		if dist == 0:
			for i, l2 in enumerate(level2):
				I_B, term, J, L, S, level, skipped_level = load_levels_lines(B_string, 1, dist, database, ele_B, charge_state, I_A)
				print("##################################################################################")
				print("##################################################################################")
				print(str(l2), str(term2[i]), str(S2[i]), str(L2[i]), str(I_B-float(l2)*hc))
				tot_cs, levels_pops_detunes, levels_cs, df_double_ce  = iter_calc(l2, S2[i], L2[i], I_B-float(l2)*hc, term, J, L, S, level, skipped_level, 1, df_double_ce, ele_B, ele_A, prod_S_A, react_S_A, react_L_A, charge_state, T_B, dist)
				print("Total CS:", tot_cs)
				second_ce_levels.append(l2)
				second_ce_cs.append(tot_cs)

			file_string='results/Z'+str(ele_B.atomic_number)+'/second_ce/levels_pops_detunes'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist)
			df_double_ce.to_csv(file_string)

			second_ce_levels_css=np.array(list(zip(second_ce_levels, second_ce_cs)), dtype=float)
			np.savetxt('results/Z'+str(ele_B.atomic_number)+'/second_ce_levels_css'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist), second_ce_levels_css, delimiter=';')
		else:
			print("can't do double cec with dist !=0 ")

	##################evolve population###################
	if not sidepeaks_sim: levels_pops_detunes=levels_pops_detunes[:,[0,1]] # remove energy differences column

	evolved_level=[]
	unevolved_level=[]

	amu=codata.value('atomic mass constant')#kg
	velocity=np.sqrt((2*T_B*1.6E-19)/(m_B*amu))*10**2
	c=codata.value('speed of light in vacuum')*10**2 #cm
	flight_time=dist/velocity
	print("flight time:", flight_time)

	print(ele_sym_A,":" ,ele_A.description)
	print(ele_sym_B,":", ele_B.description)

	E_k_probs.sort()
	E_k_probs=list(E_k_probs for E_k_probs,_ in itertools.groupby(E_k_probs)) # remove possible duplicates


	if flight_time != 0:

		dt=(flight_time/time_steps)

		for step in range(0,time_steps):
			print("time step:",step)
			#times previous step by time increment
			# cs_norm_evol_previous=list(cs_norm_evol)
			levels_pops_detunes_previous=levels_pops_detunes

			#first decrease population of the upper level energies here
			for index1,level_pop_detune1 in enumerate(levels_pops_detunes): #go through uppers
				level1=level_pop_detune1[0]#upper level
				pop1=level_pop_detune1[1]#upper pop
				if sidepeaks_sim:	eloss1=level_pop_detune1[2]

				for E_k_prob in E_k_probs:
					upper_energy=E_k_prob[0]
					lower_energy=E_k_prob[1]
					this_A=E_k_prob[2]

					if round(upper_energy,2) == round(level1, 2):#find decays from this upper, round incase diff databases

						if round(lower_energy,2) in np.around(levels_pops_detunes[:,0], decimals=2): # maintains precision of NIST database for lower energy
							index3 = np.where(np.around(levels_pops_detunes[:,0], decimals=2)==round(lower_energy,2))
							lower_energy=float(levels_pops_detunes[:,0][index3][0])

						if sidepeaks_sim:
							index2=np.array(np.where(np.all(levels_pops_detunes[:,[0,2]] == np.array([lower_energy,level1-eloss1-lower_energy]),axis=1)))

							if index2.size == 1:
								levels_pops_detunes[:,1][index2]=levels_pops_detunes[:,1][index2]+pop1*(1-np.exp(-dt*this_A)) #add pop into lower level
							elif index2.size == 0:
								newrow=np.array([lower_energy, pop1*(1-np.exp(-dt*this_A)), (level1-eloss1-lower_energy)], dtype=float)
								levels_pops_detunes=np.vstack((levels_pops_detunes,newrow))

							else: print("index 2 size error", index2)

							pop1=pop1*np.exp(-dt*this_A) #decay upper level pop
							levels_pops_detunes[:,1][index1]=pop1*np.exp(-dt*this_A)

						else:
							if lower_energy in levels_pops_detunes[:,0]:#adds pop to lower level if exiting level
								index2=np.where(levels_pops_detunes[:,0]==lower_energy)
								levels_pops_detunes[:,1][index2]=levels_pops_detunes[:,1][index2]+pop1*(1-np.exp(-dt*this_A))#adds decayed pop to lower level
							else:
								newrow=np.array([lower_energy, pop1*(1-np.exp(-dt*this_A))], dtype=float)#creates new pop if new level
								levels_pops_detunes=np.vstack((levels_pops_detunes,newrow))

							# if upper_energy==39625.506:
							# 	print("upper_energy, lower_energy, pop1, pop1*np.exp(-dt*this_A), this_A")
							# 	print(upper_energy, lower_energy, pop1, pop1*np.exp(-dt*this_A), this_A)
							# 	print(step, step*dt)

							pop1=pop1*np.exp(-dt*this_A) #decay upper level pop
							levels_pops_detunes[:,1][index1]=pop1


		print(levels_pops_detunes)
		print("sum of pops:",sum(levels_pops_detunes[:,1]))
		print("final states populated")
	else:
		print("no time of flight population change")


	if sidepeaks_sim:
		levels_pops_detunes_sorted=cumsum_diff(levels_pops_detunes[:,[0,1]]) #get rid of seperate entries for displaying populations
		levels_pops_detunes_sorted=levels_pops_detunes_sorted[levels_pops_detunes_sorted[:,1].argsort()[::-1]]
	else:
		levels_pops_detunes_sorted=levels_pops_detunes
		levels_pops_detunes_sorted=levels_pops_detunes_sorted[levels_pops_detunes_sorted[:,1].argsort()[::-1]]

	levels_pops_detunes_sorted=levels_pops_detunes_sorted.astype(float)

	print("top 10 populations:",levels_pops_detunes_sorted[:10])

	np.savetxt('results/Z'+str(ele_B.atomic_number)+'/levels_pops_detunes_sorted'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist), levels_pops_detunes_sorted, delimiter=';', fmt="%.5f")

	fig, ax = plt.subplots(figsize=(6, 10))

	# if show_notrans:
	# 	unevolved_level=list(set(energy_levels_evol)-set(evolved_level))
	# 	if len(unevolved_level)>0:
	# 		plt.axvline(x=unevolved_level[1], color='m', linestyle='-', label="No transitions from level (no info)")
	# 		for ulevel in unevolved_level[2:]: plt.axvline(x=ulevel, color='m', linestyle='-')

	if len(skipped_level)>0:
		ax.plot( [0,1] ,[skipped_level[0],skipped_level[0]], color='g', linestyle='--', label="Skipped initial population (unknown spin)")
		for slevel in skipped_level[1:]:
			ax.plot( [0,1] ,[slevel,slevel], color='g', linestyle='--')

	ax.plot(levels_pops_detunes_initial[:,1],levels_pops_detunes_initial[:,0],'ro',label="Initial population")

	for pt in levels_pops_detunes_initial[:,[0,1]]:
		ax.plot( [0,pt[1]],[pt[0],pt[0]], color='r')

	if flight_time != 0:

		ax.plot(levels_pops_detunes_sorted[:,1],levels_pops_detunes_sorted[:,0],'bs',label="Final population")
#
		for pt in levels_pops_detunes_sorted[:,[0,1]]:
			ax.plot([0,pt[1]] ,[pt[0],pt[0]], color='b')

	ax.set_ylabel("Energy level (cm-1)", fontsize=14)
	ax.set_xlabel("Normalised population", fontsize=14)
	ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	ax.plot([0,max(levels_pops_detunes_sorted[:,1])],[(I_B-I_A)/hc,(I_B-I_A)/hc], color='k', linestyle='--', label="CE entry energy")

	##resonance plot##

	if 1: #plot bare cross_section
		try:
			energy_levels, bare_cross_sections=levels_cs[:,0],levels_cs[:,1]
			bare_cs_norm_before = [float(i)/sum(bare_cross_sections) for i in bare_cross_sections] # this norm is without f factor..

			bare_cs=[float(c) for c in bare_cross_sections]
			bare_el=[float(e) for e in energy_levels]
			extra_els=np.linspace(min(bare_el),max(bare_el),30)

			for el in extra_els:
				I_B_ex=I_B-float(el)*hc
				try:
					cs=cross_section(I_A, I_B_ex, m_B, T_B)
					bare_cs.append(float(cs))
					bare_el.append(float(el))
					print("el","cs",el,cs)
				except:
					print("problem with", el, "cm-1")
				# if cs < 2.0E-10: #to reduce the amount of time wasted on calculation for the plot
				# 	break

			print("bare_cs_norm_before", bare_cs_norm_before)
			print("bare_cs", bare_cs)

			#values.index(min(values))

			bare_cs_norm_after = [float(i)/sum(bare_cs) for i in bare_cs]

			# scale_ref=bare_cs_norm_before.index(max(bare_cs_norm_before))

			# scale_ratio=bare_cs_norm_after[scale_ref]/(bare_cs_norm_before[scale_ref]*levels_pops_detunes_initial[:,1][scale_ref]) #levels_pops_detunes_initial[:,1][scale_ref] for f factor norm

			# max_pop=max(levels_pops_detunes_sorted[:,1])
			max_cs_i=bare_cs_norm_before.index(max(bare_cs_norm_before))
			# cs_i=bare_cs_norm_before[0]
			scale_ratio=max(levels_pops_detunes_sorted[:,1])/bare_cs_norm_after[max_cs_i]#)*levels_pops_detunes_initial[:,1][0]


			print("scale_ratio", scale_ratio)

			# scale_ratio=float(bare_cs_norm_before[-1])/float(bare_cross_sections[-1])

			# bare_cs_norm=[float(i)*scale_ratio for i in bare_cs]

			bare_cs_norm=[scale_ratio*float(i)/sum(bare_cs) for i in bare_cs]

			bare_el, bare_cs_norm = (list(t) for t in zip(*sorted(zip(bare_el, bare_cs_norm))))
			ax.plot(bare_cs_norm,bare_el, '--', label="Bare cross-section")

		except Exception as exception:
			print(exception)

	ax.legend(numpoints=1, loc="upper right")

	mpld3.save_html(fig=fig, fileobj='results/Z'+str(ele_B.atomic_number)+'/fig_'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist)+'.html')
	fig.savefig('results/Z'+str(ele_B.atomic_number)+'/fig_'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist)+'.pdf')
	plt.show()

	if sidepeaks_sim:
		plt.show()

	print("finished popsim")

	# change energy of beam by the energy loss of previous transitions (make sure includes chains)
	# and then dopplershift transition along with a scan range.. make last column beam energy instead?

	if sidepeaks_sim:
		try:
			E_line_rest=sidepeaks_transition[1]-sidepeaks_transition[0]
			transition_detuning_pops_f = open('results/Z'+str(ele_B.atomic_number)+'/transition_detuning_pops'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist)+"_"+str(E_line_rest),'rb')
			transition_detuning_pops=np.loadtxt(transition_detuning_pops_f, delimiter=';',dtype=float)
			# print(transition_detuning_pops)
		except Exception as exception:
			print(exception)
			#work out doppler shift of all transitions from lower levels,
			#keep/plot those with range matching near what we want and keep track of seperate transitions

			cm_to_MHz=29.9792458*10**3 #1 cm^-1=29.9792458 GHz in vacuum

			E_line_rest=sidepeaks_transition[1]-sidepeaks_transition[0] #transition cm-^1

			velocity=np.sqrt((2*(T_B)*1.6E-19)/(m_B*amu)) #shifts energy of beam in eV then converts to velocity
			beta=velocity/(c*10**(-2))
			if sidepeaks_collinear: E_line=E_line_rest*np.sqrt((1+beta)/(1-beta))
			else: E_line=E_line_rest*np.sqrt((1-beta)/(1+beta))

			only_same_transition=True

			transition_detuning_pops=[]

			for index1, level_pop_detune1 in enumerate(levels_pops_detunes): #where equals lower level..
				level1=level_pop_detune1[0]
				pop1=level_pop_detune1[1]
				eloss1=level_pop_detune1[2]

				if only_same_transition and (level1 == sidepeaks_transition[0]):

					velocity=np.sqrt((2*(T_B-eloss1*hc)*1.6E-19)/(m_B*amu)) #shifts energy of beam in eV then converts to velocity #BEAM ENERGY IN JOULES NOT eV!!
					beta=velocity/(c*10**(-2))

					if sidepeaks_collinear:
						E_line_detune=E_line_rest*np.sqrt((1+beta)/(1-beta))
					else:
						E_line_detune=E_line_rest*np.sqrt((1-beta)/(1+beta))

					detuning = (float(E_line_detune)-float(E_line))*cm_to_MHz #detuning from transition in MHz

					if abs(detuning) < 10000: #10 GHz scan range
						print([E_line,detuning,pop1])
						transition_detuning_pops.append([E_line,detuning,pop1])
					else:
						# print("detuning",abs(detuning), eloss1*hc ,"eV, too big?")
						pass

				else:

					for E_k_prob in (E_k_prob for E_k_prob in E_k_probs if round(E_k_prob[1],2) == round(level1,2)):#where lower energy is transition within 2 d.p. Finds upper levels that transition into this lower.
						this_transition=E_k_prob[0]-E_k_prob[1] #upper-lower

						velocity=np.sqrt((2*(T_B-eloss1*hc)*1.6E-19)/(m_B*amu)) #shifts energy of beam in eV then converts to velocity #BEAM ENERGY IN JOULES NOT eV!!
						beta=velocity/(c*10**(-2))

						if sidepeaks_collinear:
							this_transition_detune=this_transition*np.sqrt((1+beta)/(1-beta))
						else:
							this_transition_detune=this_transition*np.sqrt((1-beta)/(1+beta))

						detuning = (float(this_transition_detune)-float(E_line))*cm_to_MHz #detuning from transition in MHz

						if abs(detuning) < 10000: #10 GHz scan range
							print([this_transition,detuning,pop1])
							transition_detuning_pops.append([this_transition,detuning,pop1])

						else:
							# print("detuning",abs(detuning), eloss1*hc ,"eV, too big?")
							pass

			transition_detuning_pops=np.array(transition_detuning_pops)

			transition_detuning_pops[:,2]=transition_detuning_pops[:,2] / transition_detuning_pops[:,2].max(axis=0)

			np.savetxt('results/Z'+str(ele_B.atomic_number)+'/transition_detuning_pops'+ele_sym_B+"I"*charge_state+"_"+ele_sym_A+"_"+str(T_B)+"_"+str(dist)+"_"+str(E_line_rest), transition_detuning_pops, delimiter=';', fmt='%1.3f')


		print(transition_detuning_pops)
		transition_detuning_pops=transition_detuning_pops[transition_detuning_pops[:,0].argsort()] #sort into transitions

		transitions_legend={}

		for transition_detuning_pop in transition_detuning_pops:
			key=transition_detuning_pop[0]
			if key in transitions_legend:
				transitions_legend[key]=np.vstack((transitions_legend[key],np.array([transition_detuning_pop[1],transition_detuning_pop[2]]))) # updates if exists, else adds
			else:
				transitions_legend[key]=np.array([[transition_detuning_pop[1],transition_detuning_pop[2]]])

		print(transitions_legend)

		colors = cm.rainbow(np.linspace(0, 1, len(transitions_legend)))

		fig, ax = plt.subplots(figsize=(7, 5))
		axes = [ax, ax.twiny()]#, ax.twinx()

		for col_i, key in enumerate(transitions_legend):

			axes[0].plot(transitions_legend[key][:,0], transitions_legend[key][:,1], 'ro', label=str(round(key,2))+"transition", c=colors[col_i])

			for pt in transitions_legend[key][:,[0,1]]:
				# plot (x,y) pairs.
				# vertical line: 2 x,y pairs: (a,0) and (a,b)
				axes[0].plot([pt[0],pt[0]], [0,pt[1]] , color=colors[col_i])


		binning_MHz=10
		#binning sidepeak data
		no_bins=int(transition_detuning_pops[:,1].max()-transition_detuning_pops[:,1].min())/binning_MHz
		binned_detunes=scipy.stats.binned_statistic(transition_detuning_pops[:,1],transition_detuning_pops[:,2], bins=no_bins, statistic="sum")
		binned_detunes_bins=[(a + b) / 2 for a, b in zip(binned_detunes.bin_edges[0:], binned_detunes.bin_edges[1:])]

		detune_counts=[]
		for detune_i, detune in enumerate(transition_detuning_pops[:,1]):
			detune_counts.append([detune]*int(transition_detuning_pops[:,2][detune_i]*1000))

		# detune_counts=list(np.array(detune_counts).flat)
		# print("detune_counts", detune_counts)
		detune_counts=[item for sublist in detune_counts for item in sublist]

		#binning sidepeak data
		axes[0].bar(binned_detunes_bins, binned_detunes.statistic, label=str(binning_MHz)+" MHz binning", align='center', alpha=0.4, width=(binned_detunes_bins[1]-binned_detunes_bins[0])) # A bar chart

		# X_plot = np.linspace(min(transition_detuning_pops[:,1])-50, max(transition_detuning_pops[:,1])+50, 1000)
		# kde = KernelDensity(kernel='gaussian', bandwidth=4).fit(np.array(detune_counts).reshape((-1, 1)))
		# log_dens = kde.score_samples(np.array(X_plot).reshape((-1, 1)))
		# plt.fill(X_plot, np.exp(log_dens)*1000.0/(35*1.75), fc='#AAAAFF')

		# sidepeak_model=ConstantModel()
		# for detune_i, detune in enumerate(transition_detuning_pops[:,1]):
		# 	sidepeak_model=sidepeak_model + LorentzianModel() # amplitude=np.array([transition_detuning_pops[:,2][detune_i]], sigma=[10.0], center=[detune]
		#
		# print(sidepeak_model.param_names)
		#
		# print("X_plot, sidepeak_model(X_plot)", X_plot, sidepeak_model(X_plot))
		# plt.plot(X_plot, sidepeak_model(X_plot))

		axes[0].set_xlabel("Detuning (MHz)", fontsize=20)
		axes[1].set_xlabel('Detuning (eV)', fontsize=20)
		plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		plt.ylabel("Norm. Population", fontsize=20)
		plt.legend(numpoints=1, loc=2)
		plt.xlim((-200, +200))

		# ax1Ticks = axes[0].get_xticks()
		# ax2Ticks = ax1Ticks
		#
		# def tick_function(X):
		# 	V = X + 49
		# 	return [str(z) for z in V]
		#
		# axes[1].set_xticks(ax2Ticks)
		# axes[1].set_xbound(axes[0].get_xbound())
		# axes[1].set_xticklabels(tick_function(ax2Ticks))

		ax = plt.gca()
		ax.get_xaxis().get_major_formatter().set_scientific(False)

		plt.show()

with open(sim_list,'rb') as csv_file:

	iter_list = np.genfromtxt(csv_file,delimiter = '	',dtype=str,skip_header=1,autostrip=1)
	sim_no=[int(x) for x in iter_list.T[0]]
	ele_sym_Bs=[x for x in iter_list.T[1]]
	ele_sym_As=[ast.literal_eval(x) for x in iter_list.T[2]]
	sidepeaks_transitions=[ast.literal_eval(x) for x in iter_list.T[3]]
	sidepeaks_collinears=[ast.literal_eval(x) for x in iter_list.T[4]]
	dists=[float(x) for x in iter_list.T[5]]
	T_Bs=[float(x) for x in iter_list.T[6]]
	time_stepss=[int(x) for x in iter_list.T[7]]
	skips=[bool(int(x)) for x in iter_list.T[8]]
	databases=[str(x) for x in iter_list.T[9]]
	charge_states=[int(x) for x in iter_list.T[10]]
	double_ces=[bool(int(x)) for x in iter_list.T[11]]
	csv_file.close()


for sym_B_i, ele_sym_B in enumerate(ele_sym_Bs):
	for sym_A_i, ele_sym_A in enumerate(ele_sym_As[sym_B_i]):

		skip=skips[sym_B_i]
		database=databases[sym_B_i]
		sidepeaks_transition=sidepeaks_transitions[sym_B_i]
		sidepeaks_collinear=sidepeaks_collinears[sym_B_i]
		dist=dists[sym_B_i]
		T_B=T_Bs[sym_B_i]
		time_steps=time_stepss[sym_B_i]
		charge_state=charge_states[sym_B_i]
		double_ce=double_ces[sym_B_i]

		print("tion:",ele_sym_A,ele_sym_B,T_B,dist,sidepeaks_transition,sidepeaks_collinear,time_steps, charge_state)

		if not skip:
			try:
				popsim(ele_sym_A,ele_sym_B,T_B,dist,sidepeaks_transition,sidepeaks_collinear,time_steps,database, charge_state, double_ce)
			except Exception as e:
				print(e)
				print("PROBLEM with element", ele_sym_B)
		else:
			print("simulation skipped")

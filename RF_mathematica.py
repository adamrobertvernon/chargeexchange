from scipy.constants import codata
import scipy as sp
from timeout import *
import numpy as np
from subprocess import *

def runMath(parameter):
	command='/usr/bin/runMath'
	output = check_output([command,parameter])
	output=output.decode('ascii')
	output=output.replace("*^","E")
	try:
		output=output.split('\n')#splits lines into items
		output=float(output[-2])#gets last line returned from Mathematica even if there are leading error messages
		return output
	except:
		print(output)

@timeout(30)
def cross_section(I_A,I_B,m_B,T_B):
	#a_0=codata.value('Bohr radius')#m
	a_0=5.291E-9 #cm
	hbar=6.582E-16 #eV s
	m_p=sp.constants.proton_mass
	m_n=sp.constants.neutron_mass
	m_e=sp.constants.electron_mass
	pi=sp.constants.pi
	e=sp.constants.elementary_charge
	k=sp.constants.Boltzmann
	amu=codata.value('atomic mass constant')#kg

	Ibar=(I_A+I_B)/2 #eV
	if Ibar<0: return 0.0 #it would be unbound/continuum
	DeltaE=abs(I_A-I_B) #eV
	omega=DeltaE/hbar #/s
	gamma=np.sqrt(Ibar/13.6)
	f=1

	v=np.sqrt((2*T_B*1.6E-19)/(m_B*amu))*10**2

	a_0=format(a_0, '.12f')
	hbar=format(hbar, '.19f')
	omega=format(omega, '.2f')
	gamma=format(gamma, '.19f')
	Ibar=format(Ibar, '.19f')
	v=format(v, '.19f')

	s1=("(Sin[Sqrt[(2*Pi)/("+str(gamma)+"*"+str(a_0)+")]*(2*"
	""+str(Ibar)+")/("+str(hbar)+"*"+str(v)+")"
	"*b^(3/2)*(1 + "+str(a_0)+"/("+str(gamma)+"*b))*"
	"Exp[(-"+str(gamma)+"*b)/"+str(a_0)+"]]^2)*"
	"(Sech[("+str(omega)+"/"+str(v)+")*Sqrt[("+str(a_0)+"*Pi*b)/(2*"+str(gamma)+")]]^2)*2*Pi*b")

	#if calculation of extremely small cross-sections is needed then MinRecursion can be increased and timeout set to 0
	s2="Sqrt["+str(f)+"*NIntegrate["+s1+", {b, 0, Infinity}, MaxRecursion -> 200, MinRecursion -> 10]]"

	print(gamma,Ibar,omega,v)

	result = runMath(s2)**2

	if result != None:
		return result
	else:
		return 0.0

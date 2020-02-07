# 2020-02-04,1.15pm: Started synapse_groups branch. 
# Objective: refine the two-synapse group architecture
# 2020-02-07: add IP3 and Ce equations, even though the code is not certified to be correct. 

#  SIMPLIFY THE CODE TO ONLY KEEP THE Calcium and look at diffusion term only. 

from brian2 import *
import utils as u
#codegen.target = 'cython'

### General parameters
# if dt=1*ms, I get Nan
# VVT is very non-smooth as a function of time
# Following Maurizio's derivation, implement the solution to the cubic equation in Cytosolic Calcium concentration.
# Otherwise very similar to voltage_cubic.py

# Initial conditions should be self-consistent in some way. 
# What is V0, C0, etc at t=0?

duration = 0.7*ms          # Total simulation time
sim_dt = 0.1*ms            # Integrator/sampling step

# Model definition
defaultclock.dt = sim_dt     # Set the integration time

# Extracellular calcium (for Goldman formula). 
# What is the location? Where was Goldman's formula applied? Across the extracelluar membrane (radius R), 
# or between the cytosol and the ER? 
Ce  = 1.e-4 * umole / cm**3
print("Ce= ", Ce)
Ce  = 1800 * uM   # too large? 
print("Ce= ", Ce)
# From Oschmann Master thesis, page 53. 
Crest = 0.073 * uM  # M is molar, u is micro
print("Crest= ", Crest)

P_Ca = 4.46e-13 * cm / second   # From "Calcium Permeability in Large Unilamellar Vesicles Prepared from Bovine Lens Cortical Lipids"
V_T = 26 * mvolt
Vrest  = -80 * mvolt
F = 96485.3329 * amp * second / mole  #  is Coulomb
k_B = 1.38064852e-23 * meter**2 * kilogram / second**2 / kelvin
N_A = 6.02214076e23 / mole
e = 1.602176634e-19 * amp * second
Cm = 1 * uF / meter**2 #: 1  #Capacitance per meter**2
D_C  = 5.3e-6 * cm**2 / second # from "Free diffusion coefficient of ionic calcium in cytoplasm", Donahue et al (1987)
D_CE = 5.3e-6 * cm**2 / second # (a guess, ER)
D_I  = 5.3e-6 * cm**2 / second # (a guess, IP3)
Rgas = 8.31 * joule / kelvin / mole

#----------------------------------------------------------------------
astro_eqs = '''
s = R + r          : meter
r                  : meter
dR2 = R**2 - r**2  : meter**2
DR2                : meter**2
L                  : meter
R                  : meter  # Outer radius
coupling_C         : mole/second/meter**3

####  CALCIUM

# has large impact on diffusion. Creates instability
coupling_electro   : mole / second / meter**3
#coupling_electro = 1  : 1 

# Tot_C is the total calcium (computed in synapses2)
Tot_C : mole/meter**3

V = Vrest + (F*s/Cm) * (Tot_C-Crest) : volt
VVT = -2.*V/V_T : 1  # value of 0, which leads to a zero denominator

# Has minimal impact on Calcium
electro_diffusion = -P_Ca * V / (R * V_T) * (Ce*exp(-2*V/V_T) - Tot_C) / (exp(-2*V/V_T) - 1) : mole / second / meter**3

A1 : meter**4 / mole 
B1 : meter 
C1 : mole / meter**2 

CC0 : mole / meter**3
V0 : volt
nb_connections : 1  # number of synaptic connections

# Calcium diffusion
# Each term will handle half a compartment 

dC/dt = 1*coupling_C + 1.*coupling_electro + 0.*electro_diffusion        : mole / meter**3

####  ENDOPLASMIC RETICULUM
Tot_CE                 : mole/meter**3
coupling_CE            : mole/second/meter**3
dCE/dt = 1*coupling_CE : mole/meter**3

####  IP3
Tot_I                  : mole/meter**3
coupling_I             : mole/second/meter**3
dI/dt = 1*coupling_CE  : mole/meter**3

### Open Channels
dh/dt = OmegaH * (Hinf - h) : 1


################3
###  CURRENTS

Jr     =
J1     = 
Jp     = 
minf   = 
ninf   = 
hinf   = 
OmegaH = 
Jbeta  = 
Jdelta = 
J3K    = 
J5P    = 



#IP3 dynamics
Jdelta=O_delta * ( 1 - ( IP3 /(kappa_delta + IP3) ) * ( CaCYT**2/(CaCYT**2+K_delta**2) )  ) : mole/second
J3k=O_3K * ( ( CaCYT**4/(CaCYT**4+K_D**4)) * (IP3/(K_3K+IP3) ) ) : mole/second
J5p=O_5P * IP3 /(IP3 + K_5P) : mole/second

# Calcium core
minf=IP3/(IP3+d_1) : 1
ninf=CaCYT/(CaCYT+d_5) : 1
Jcicr=Omega_C*(minf*ninf*h)**3*(CaER-CaCYT) : mmolar/second
Jleak=Omega_L*(CaER-CaCYT) : mmolar/second
Jserca=O_P*CaCYT**2/(K_P**2+CaCYT**2) : mmolar/second
'''

N_astro = 3 # Total number of astrocytes1 in the network (each compartment broken into two
N_astro = 2*N_astro # Total number of astrocytes1 in the network (each compartment broken into two
astrocytes1 = NeuronGroup(N_astro, astro_eqs, method='euler', name='ng1', order=1)

# Initital Conditions
astrocytes1.L = 8 *  umeter  # length of a compartment
astrocytes1.R = 3 * umeter
astrocytes1.r = 2 * umeter
# Calcium in 1/2 compartment should be the value for the full compartment since the 1/2 compartment
# Each 1/2 compartment only considers the fork at one end. 
# is an artifact of the numerical scheme.
astrocytes1.C = 1.1e-4  * mole / meter**3
astrocytes1.I = 1.1e-4  * mole / meter**3
astrocytes1.CE = 1.1e-4  * mole / meter**3
#astrocytes1.C = [1.1e-4, 1.5e-4, 1.6e-4]*2 * mole / meter**3

astro_mon = StateMonitor(astrocytes1, variables=['Tot_C', 'Tot_CE', 'Tot_I', 'coupling_C', 'coupling_electro', 'electro_diffusion', 'nb_connections', 'A1', 'B1', 'C1', 'CC0', 'V0'], record=True)
b_mon = StateMonitor(astrocytes1, variables=['CC0', 'dR2', 's', 'L'], record=True)

'''
One issue is Boundary Conditions. I have not done anything to implement those. 
Also, each synapse has its own value of V0. 
'''

#----------------------------------------------------------------------
# Diffusion between astrocytes1
synapse1_eqs = '''
# Assuming the above is correct, let us figure out units. 
# using C_pre:  C * F * s / Cm = V  ==> F / Cm = V / s / C = (V / s) * m^3 / mole

####  CYTOPLASM

# A1 stable: value does not change from iteration to iteration
A1_pre = (0.5 * F * s_post)/(Cm * V_T) * (dR2_post / L_post) : meter**4 / mole (summed)
B1_pre = (1. - (s_post * F * Tot_C_post)/(2. * Cm * V_T)) : meter (summed)
C1_pre = dR2_post * Tot_C_post / L_post : mole / meter**2 (summed)

# Only computed to trick Brian into computing C0_pre by summing the values for all synapses (which will 
#   all bequal  since it depends on A1, B1, C1, for which the summation has already been executed. 
# There must be an easier method without calculating C0 and V0 in the NeuronGroup
nb_connections_pre = 1 : 1 (summed)

# Notice I am dividing by the number of connections (3 if a single fork)
CC0_pre = ((-B1_pre + sqrt(B1_pre**2 + 4*C1_pre*A1_pre)) / (2.*A1_pre)) / nb_connections_pre : mole / meter**3 (summed)
V0_pre = (Vrest + (CC0_pre - Crest) * (F * s) / Cm) / nb_connections_pre   : volt  (summed)

# Investigate why this should be a minus sign
coupling_C_post = (4*D_C / L_post**2) * (Tot_C_pre - Tot_C_post) : mole/second/meter**3 (summed)

# MAKE SURE CC0 and V0 are computed before updating the coupling parameters. 
coupling_electro_pre = (4*D_C/L_post**2/V_T) * (CC0_post + Tot_C_post) * (V0_post - V_post) : mole/second/meter**3 (summed) 

####  ENDOPLASMIC RETICULUM
coupling_CE_post = (4*D_CE / L_post**2) * (Tot_CE_pre - Tot_CE_post) : mole/second/meter**3 (summed)

####  IP3
coupling_I_post = (4*D_I / L_post**2) * (Tot_I_pre - Tot_I_post) : mole/second/meter**3 (summed)
'''
#----------------------------------------------------------------------

# TEMPORARY
synapse2_eqs = '''
    # C is updated in each compartment of the astrocyte and combined here
	Tot_C_syn   = C_pre  + C_post    : mole / meter**3
	Tot_C_pre   = Tot_C_syn       : mole / meter**3 (summed)

	# update Calcium in compartments
	# This works because synapses are bidirectional i --> j and j --> i

####  ENDOPLASMIC RETICULUM
	Tot_CE_pre  = Tot_CE_syn      : mole / meter**3 (summed)
	Tot_CE_syn  = CE_pre + CE_post   : mole / meter**3

####  IP3
	Tot_I_syn = I_pre + I_post : mole / meter**3
	Tot_I_pre = Tot_I_syn      : mole / meter**3 (summed)
'''

synapses1 = Synapses(astrocytes1, astrocytes1, model=synapse1_eqs, method='euler', order=0, name='sg1')
synapses2 = Synapses(astrocytes1, astrocytes1, model=synapse2_eqs, method='euler', order=5, name='sg2')

# Connections count from 0
#  ---*---.---*----
#         |___*____

# Single fork (N_astro = 6)
# Six Compartments: 0,1,...,4,5
for i in range(N_astro):
    synapses1.connect(i=i,j=i)
pairs = [(3,1),(3,2),(1,2)]
for pair in pairs:
	synapses1.connect(i=pair[0], j=pair[1])
	synapses1.connect(i=pair[1], j=pair[0])

pairs = [(0,1),(2,4),(3,5)]
for pair in pairs:
	synapses2.connect(i=pair[0], j=pair[1])
	synapses2.connect(i=pair[1], j=pair[0])

s = synapses1.summed_updaters
s['nb_connections_pre'].order = 0
s['A1_pre'].order = 2
s['B1_pre'].order = 2
s['C1_pre'].order = 2
s['CC0_pre'].order = 4
s['V0_pre'].order = 5
s['coupling_C_post'].order = 6
s['coupling_electro_pre'].order = 7

astrocytes1.state_updater.order = 10

s = synapses2.summed_updaters
s['Tot_C_pre'].order = 13


for k,v in s.items():
	print("updater[%s]= " % k, v)


# Must call StateMonitor AFTER synaptic connections are established
syn1_mon = StateMonitor(synapses1, variables=[], record=True, name='syn1_mon')
syn2_mon = StateMonitor(synapses2, variables=[], record=True, name='syn2_mon')

groups = [synapses1, synapses2]  #  must have i,j connections
u.connectionMatrices(groups)


#-------------------------------------------
# Initial Conditions (call after connect())
# Temporary. Should be computed in the synapse
#synapses1.P                  = 1.1e-2 

groups = [astrocytes1, synapses1, synapses2]
u.lowLevelInfo(groups)

# Run Simulation
print(scheduling_summary())
run(duration)
#run(duration, report='text')

u.printData(astro_mon, b_mon, syn1_mon, syn2_mon)
quit()
#u.plots(astro_mon)
#----------------------------------------------------------------------


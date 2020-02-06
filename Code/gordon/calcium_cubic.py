# 2020-02-04,1.15pm: Started synapse_groups branch. 
# Objective: refine the two-synapse group architecture

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

duration = 0.3*ms          # Total simulation time
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
D_C = 5.3e-6 * cm**2 / second # from "Free diffusion coefficient of ionic calcium in cytoplasm", Donahue et al (1987)
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

# has large impact on diffusion. Creates instability
coupling_electro   : mole / second / meter**3
#coupling_electro = 1  : 1 

# Tot_C is the total calcium (computed in synapses2)
Tot_C : mole/meter**3

V = Vrest + (F*s/Cm) * (Tot_C-Crest) : volt
VVT = -2.*V/V_T : 1  # value of 0, which leads to a zero denominator

# Has minimal impact on Calcium
electro_diffusion = -P_Ca * V / (R * V_T) * (Ce*exp(-2*V/V_T) - Tot_C) / (exp(-2*V/V_T) - 1) : mole / second / meter**3
#electro_diffusion :  mole / second / meter**3  # computed in synapse
#electro_diffusion = 1 : 1

A1 : meter**4 / mole 
B1 : meter 
C1 : mole / meter**2 

CC0 : mole / meter**3
V0 : volt
A000_nb_connections : 1  # number of synaptic connections
nb_connections2 : 1  # number of synaptic connections in synapses2

# Calcium diffusion
# Each term will be divided into two and recombined in the second synapse group into Tot_C

dC/dt = coupling_C + 0.*coupling_electro + 0.*electro_diffusion        : mole / meter**3


'''

N_astro = 3 # Total number of astrocytes1 in the network (each compartment broken into two
N_astro = 2*N_astro # Total number of astrocytes1 in the network (each compartment broken into two
astrocytes1 = NeuronGroup(N_astro, astro_eqs, method='euler', order=0, name='ng1')

# Initital Conditions
astrocytes1.L = 8 *  umeter  # length of a compartment
astrocytes1.R = 3 * umeter
astrocytes1.r = 2 * umeter
astrocytes1.C = 1.1e-4  * mole / meter**3
#astrocytes1.C = [1.1e-4, 1.5e-4, 1.6e-4]*2 * mole / meter**3

astro_mon = StateMonitor(astrocytes1, variables=['Tot_C', 'coupling_C', 'coupling_electro', 'electro_diffusion', 'nb_connections2', 'A000_nb_connections', 'A1', 'B1', 'C1', 'CC0', 'V0'], record=True)
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

dP/dt = 1*Hz : 1

# Solution to quadratic C0 equation (calcium concentration)
# A1 C0^2 + B1 C - C1 = 0
# C0 = (-B1 \pm \sqrt{B1^2 + 4*C1*A1) / (2*A1)

# A1 stable: value does not change from iteration to iteration
A1_pre = (0.5 * F * s_post)/(Cm * V_T) * (dR2_post / L_post) : meter**4 / mole (summed)
B1_pre = (1. - (s_post * F * Tot_C_post)/(2. * Cm * V_T)) : meter (summed)
C1_pre = dR2_post * Tot_C_post / L_post : mole / meter**2 (summed)

# Only computed to trick Brian into computing C0_pre by summing the values for all synapses (which will 
#   all bequal  since it depends on A1, B1, C1, for which the summation has already been executed. 
# There must be an easier method without calculating C0 and V0 in the NeuronGroup
A000_nb_connections_pre = 1 : 1 (summed)

# Figure out which branch to use: +sqrt or -sqrt? 
# or -sqrt... But C0 > 0 (is C1*A1 always greater than 0?)

# Notice I am dividing by the number of connections (3 if a single fork)
CC0_pre = ((-B1_pre + sqrt(B1_pre**2 + 4*C1_pre*A1_pre)) / (2.*A1_pre)) / A000_nb_connections_pre : mole / meter**3 (summed)
V0_pre = (Vrest + (CC0_pre - Crest) * (F * s) / Cm) / A000_nb_connections_pre   : volt  (summed)

# Investigate why this should be a minus sign
coupling_C_post = (4*D_C / L_post**2) * (Tot_C_pre - Tot_C_post) : mole/second/meter**3 (summed)

# MAKE SURE CC0 and V0 are computed before updating the coupling parameters. 
coupling_electro_pre = (4*D_C/L_post**2/V_T) * (CC0_post + Tot_C_post) * (V0_post - V_post) : mole/second/meter**3 (summed) # diverge
'''
#----------------------------------------------------------------------

# TEMPORARY
synapse2_eqs = '''
    # C is updated in each compartment of the astrocyte and combined here
	Tot_C_syn = C_pre + C_post : mole / meter**3

	# update Calcium in compartments
	# This works because synapses are bidirectional i --> j and j --> i
    nb_connections2_pre = 1 : 1 (summed)

	# How many connection? 
	Tot_C_pre   = Tot_C_syn       : mole / meter**3 (summed)
	#Tot_C_post  = Tot_C_syn      : mole / meter**3 (summed)
'''

synapses1 = Synapses(astrocytes1, astrocytes1, model=synapse1_eqs, method='euler', order=2, name='sg1')
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

# Must call StateMonitor AFTER synaptic connections are established
syn1_mon = StateMonitor(synapses1, variables=['P'], record=True, name='syn1_mon')
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


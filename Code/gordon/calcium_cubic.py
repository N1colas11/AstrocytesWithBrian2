
#  SIMPLIFY THE CODE TO ONLY KEEP THE Calcium and look at diffusion term only. 

from brian2 import *
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

astro_eqs = '''
s = R + r          : meter
r                  : meter
dR2 = R**2 - r**2  : meter**2
DR2                : meter**2
L                  : meter
R                  : meter  # Outer radius
coupling_C         : mole / second / meter**3

# has large impact on diffusion. Creates instability
coupling_electro   : mole / second / meter**3

V = Vrest + (F*s/Cm) * (C-Crest) : volt
VVT = -2.*V/V_T : 1  # value of 0, which leads to a zero denominator

# Has minimal impact on Calcium
electro_diffusion = -P_Ca * V / (R * V_T) * (Ce*exp(-2*V/V_T) - C) / (exp(-2*V/V_T) - 1) : mole / second / meter**3

A1 : meter**4 / mole 
B1 : meter 
C1 : mole / meter**2 

CC0 : mole / meter**3
V0 : volt
A000_nb_connections : 1  # number of synaptic connections

# Calcium diffusion
# C: [[ 1.10000000e-07 -3.80029781e+07             nan
dC/dt = coupling_C + 1.*coupling_electro + 0.*electro_diffusion        : mole / meter**3
'''

N_astro = 3 # Total number of astrocytes1 in the network
astrocytes1 = NeuronGroup(N_astro, astro_eqs, method='euler', order=0, name='ng1')

# Initital Conditions
astrocytes1.L = 8 *  umeter  # length of a compartment
astrocytes1.R = 3 * umeter
astrocytes1.r = 2 * umeter
astrocytes1.C = [1.1e-4, 1.5e-4, 1.6e-4] * mole / meter**3
astrocytes1.C = [1.1e-4, 1.1e-4, 1.6e-4] * mole / meter**3

astro_mon = StateMonitor(astrocytes1, variables=['C','coupling_C', 'coupling_electro', 'electro_diffusion', 'A000_nb_connections', 'A1', 'B1', 'C1', 'CC0', 'V0'], record=True)
b_mon = StateMonitor(astrocytes1, variables=['VVT','V', 'CC0', 'dR2', 's', 'L'], record=True)

# Diffusion between astrocytes1
synapse1_eqs = '''
# Assuming the above is correct, let us figure out units. 
# using C_pre:  C * F * s / Cm = V  ==> F / Cm = V / s / C = (V / s) * m^3 / mole

# Solution to quadratic C0 equation (calcium concentration)
# A1 C0^2 + B1 C - C1 = 0
# C0 = (-B1 \pm \sqrt{B1^2 + 4*C1*A1) / (2*A1)

A1_pre = (0.5 * F * s_post)/(Cm * V_T) * (dR2_post / L_post) : meter**4 / mole (summed)
B1_pre = (1. - (s_post * F * C_post)/(2. * Cm * V_T)) : meter (summed)
C1_pre = dR2_post * C_post / L_post : mole / meter**2 (summed)

# Only computed to trick Brian into computing C0_pre by summing the values for all synapses (which will 
#   all bequal  since it depends on A1, B1, C1, for which the summation has already been executed. 
# There must be an easier method without calculating C0 and V0 in the NeuronGroup
A000_nb_connections_pre = 1 : 1 (summed)

# Figure out which branch to use: +sqrt or -sqrt? 
# or -sqrt... But C0 > 0 (is C1*A1 always greater than 0?)

# worked. Nb_connections was nonzero, equal to 3, as expected. When I remove it, CC0 is three times higher.
CC0_pre = ((-B1_pre + sqrt(B1_pre**2 + 4*C1_pre*A1_pre)) / (2.*A1_pre)) / A000_nb_connections_pre : mole / meter**3 (summed)
V0_pre = (Vrest + (CC0_pre - Crest) * (F * s) / Cm) / A000_nb_connections_pre   : volt  (summed)

# Investigate why this should be a minus sign
coupling_C_post = (4*D_C / L_post**2) * (C_pre - C_post) : mole / second / meter**3 (summed)

# MAKE SURE CC0 and V0 are computed before upding the coupling parameters. 
#coupling_electro_pre = (4*D_C/L_post**2/V_T) * (CC0_post + C_post) * (V0_post - V_post) : mole/second/meter**3 (summed) # diverge

# Same divergence whether _pre or _post. WHY?   (CREATES INSTABILITY after 1-2 times steps)
coupling_electro_pre = (4*D_C/L_post**2/V_T) * ( (CC0_post + C_post) * (V0_post - V_post)  - 
   (C_pre + C_pre) * (V0_pre - V_pre) )  : mole/second/meter**3 (summed)

#coupling_electro_pre = (4*D_C/L_post**2/V_T) * (C0_post + C_post) * (V0_post - V_post) : mole/second/meter**3 (summed)
'''


'''
Neuron1 <-- Synapse1 --> Neuron2 <-- Synapse2 --> Neuron1 

Mail to Maurizio, 2020-01-25,3.37pm
I am thinking that I need connections similar to:

Neuron1 <-- Synapse1 --> Neuron2 <-- Synapse2 --> Neuron1

connect()  (all to all) feels complicated. Although all to all might work with 3 neurons, it will not work with more than 3.
But this structure is required to control the computation of intermediate variables.

Neuron1 would advance C, Ce and IP3, also V, which only depends on C ( and constants Crest and Vrest)

Neuron2 would compute V0 and C0 (which have to be callable from Neuron1). Therefore, Neuron1[i] must be the same
compartment as Neuron2[i] (I do not know whether the notation is correct). Finally, Neuron2 would advance C, Ce, IP3

Synapse1 would compute A, B, C, necessary to solve the quadratic equation allowing the computation of C0 and V0 in Neuron2.
A C0**2 + B * C0 + C = 0   ,   solve for C0 (quadratic equation in C0)

Synapse2 would compute the diffusion and electro-coupling terms, which depend on  V0 and C0 (computed in Neuron2). This term also depends on C and V, which are available from Neuron1. So Synapse 2 must access variables in Neuron1 and Neuron2). I wonder what that does to efficiency.

What do you think?

Consider doing:  (use an existing structure to set the connection)
A = randint(2,size=(4,4))
rows, cols = nonzero(A)
S.connect(rows, cols)  # or cols, rows...

-------

An alternative approach would be to have a single NeuronGroup and a single SynapseGroup. In that case, each group is divided into two parts .We create a schedule executed in the following order:

Synapse1 (part1), Neuron1 (part1), Synapse1 (part2), Neuron2 (part2).

Would it be possible to add a switch controlled by a global variable. If the switch is True, part 1 is executed. If switch is False, part 2 is executed. I must be able to control the boolean value of the switch. There are two switches actually. One for the SynapseGroup, and one for the NeuronGroup.

What do you think? Is it possible?
'''

synapses1 = Synapses(astrocytes1, astrocytes1, model=synapse1_eqs, method='euler', order=2, name='sg1')
synapses1.connect()
matrix = np.zeros([len(astrocytes1), len(astrocytes1)], dtype=bool)
matrix[synapses1.i[:], synapses1.j[:]] = True
print("Connection matrix:")
print(matrix)

#syn_mon = StateMonitor(synapses1, variables=['A000_nb_connections'], record=True)

#for o in net.objects: 
	#print("object: ", o)

#for s in net.get_states():
    #print("state: ", s)

#synapses1.variables['A000_nb_connections_pre'].name = 'A000_nb_connections_pre'
#print("synapse1.variables: ", synapses1.variables['A000_nb_connections_pre'].name)
#quit()

for k,v in synapses1.variables.items():
	#print("k: ", k)
	#print("v: ", v)
	#print("..dir(k)= ", dir(k))
	#print("..dir(v)= ", dir(k))
	pass


print("..dir(synapses1)= ", dir(synapses1))
print("..dir(synapses1.subexpression_updater)= ", dir(synapses1.subexpression_updater))
print("..dir(synapses1.summed_updaters)= ", dir(synapses1.summed_updaters))
print()

# Run Simulation
print(scheduling_summary())
run(duration, report='text')

print("C= ", astro_mon.C)
print("V= ", b_mon.V)  
#print("CC0= ", b_mon.CC0)
#print("coupling_C= ", astro_mon.coupling_C)
#print("coupling_electro= ", astro_mon.coupling_electro)
#print("electro_diffusion= ", astro_mon.electro_diffusion)
#print("VVT=-2*V/V_T ", b_mon.VVT)
#print("V0= ", astro_mon.V0)
#print("V_T= ", V_T)
print("A000_nb_connections: ", astro_mon.A000_nb_connections)
print("A1= ", astro_mon.A1)
print("B1= ", astro_mon.B1)
print("C1= ", astro_mon.C1)
print("CC0= ", astro_mon.CC0)
print("V0= ", astro_mon.V0)
#print("s= ", b_mon.s)
print("dR2= ", b_mon.dR2)
print("L= ", b_mon.L)



################################################################################
# Analysis and plotting
################################################################################
fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, figsize=(6.26894 * 2, 6.26894 * 0.66),
                       gridspec_kw={'left': 0.1, 'bottom': 0.12})

ax1.plot(astro_mon.t/second,astro_mon.C[0]/umolar,label='C1')
ax1.plot(astro_mon.t/second,astro_mon.C[1]/umolar,label='C2')
ax1.legend()
ax1.set(xlabel='time (ms)',ylabel='Ca concentration (umolar)')

plt.savefig("plot.pdf")
#plt.show()

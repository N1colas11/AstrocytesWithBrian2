from brian2 import *
#codegen.target = 'cython'

### General parameters
duration = 0.3*ms       # Total simulation time
sim_dt = 0.1*ms            # Integrator/sampling step

# Model definition
defaultclock.dt = sim_dt     # Set the integration time

Ce  = .5 * mole / meter**3
P_Ca = 1. * cm / second  # Permeability
V_T = 25 * mvolt
Crest  = 1.e-4 * mole / meter**3
Vrest  = 20 * mvolt
F = 96485.3329 * second * amp / mole
Cm = 1 * uF / meter**2 #: 1  #Capacitance per meter**2
D_C = 1 * meter**2 / second

astro_eqs = '''
g                  : 1
#nb_connections    : 1
DR2                : meter**2
L                  : meter
R                  : meter  # Outer radius
s                  : meter
A                  : 1/volt
B                  : 1
CC                 : volt
dR2                : meter**2
coupling           : mole / second / meter**3
coupling_electro   : mole / second / meter**3

V = Vrest + (F*s/Cm) * (C-Crest) : volt


VVT = -2.*V/V_T : 1  # value of 0, which leads to a zero denominator
electro_diffusion = -P_Ca * V / (R * V_T) * (Ce*exp(-2*V/V_T) - C) / (exp(-2*V/V_T) - 1) : mole / second / meter**3
#electro_diffusion = -P_Ca * V / (R * V_T) * (Ce*exp(-2*V/V_T) - C) : mole / second / meter**3

C0 = Crest + (V0 - Vrest) * Cm / (F * s)   : mole / meter**3
# V0 is the problem
V0 = (B + sqrt(B**2 - 4.*A*CC)) / A        : volt   # at compartment edges

# Calcium diffusion
dC/dt = coupling + coupling_electro + electro_diffusion        : mole / meter**3
'''

N_astro = 3 # Total number of astrocytes in the network
astrocytes = NeuronGroup(N_astro, astro_eqs)
astrocytes.g = [10,20,40]
astrocytes.L = 1 *  umeter
astrocytes.s = 1 *  umeter
astrocytes.dR2 = 1 *  umeter**2
astrocytes.R = 2 * umeter
astrocytes.A = 1 / volt


astro_mon = StateMonitor(astrocytes, variables=['C','coupling', 'coupling_electro', 'electro_diffusion', 'A', 'V0'], record=True)
b_mon = StateMonitor(astrocytes, variables=['VVT','V'], record=True)

# Diffusion between astrocytes
astro_to_astro_eqs = '''
A_pre = dR2_post / (L_post * s_post * V_T) : 1/volt (summed)
B_pre = 0.5 * (dR2_post / (L_post * s_post)) * (Vrest + V_post)/V_T : 1 (summed)
CC_pre = (dR2_post/L_post * s_post) * ((Crest-C_post) * (F * s_post/Cm) + Vrest * (V_post/V_T-1.)) : volt (summed)

# The next two terms appear to be equal once summed. Strange.

coupling_pre = (4*D_C / L_post**2) * (C0_post - C_post) : mole / second / meter**3 (summed)

coupling_electro_pre = (4*D_C/L_post**2/V_T) * (C0_post + C_post) * (V0_post - V_post) : mole/second/meter**3 (summed)

'''

astro_to_astro = Synapses(astrocytes, astrocytes, model=astro_to_astro_eqs)
#syn_mon = StateMonitor(astro_to_astro, variables=['fff'], record=True)
astro_to_astro.connect()

# Simulation run
run(duration, report='text')

print("C= ", astro_mon.C)
#print("ddd= ", astro_mon.ddd)
#print("eee= ", astro_mon.eee)
print("coupling= ", astro_mon.coupling)
print("coupling_electro= ", astro_mon.coupling_electro)
print("electro_diffusion= ", astro_mon.electro_diffusion)
#print("nb_connections= ", astro_mon.nb_connections)
#print("fff= ", syn_mon.fff)
print("A= ", astro_mon.A)
print("VVT=-2*V/V_T ", b_mon.VVT)
print("V0= ", astro_mon.V0)
print("V= ", b_mon.V)
print("V_T= ", V_T)

from brian2 import *
#codegen.target = 'cython'

### General parameters
# if dt=1*ms, I get Nan
# VVT is very non-smooth as a function of time
duration = 30*ms          # Total simulation time
sim_dt = .1*ms            # Integrator/sampling step

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
#nb_connections    : 1
DR2                : meter**2
L                  : meter
R                  : meter  # Outer radius
A                  : 1/volt
B                  : 1
CC                 : volt
coupling           : mole / second / meter**3
coupling_electro   : mole / second / meter**3

V = Vrest + (F*s/Cm) * (C-Crest) : volt


VVT = -2.*V/V_T : 1  # value of 0, which leads to a zero denominator
electro_diffusion = -P_Ca * V / (R * V_T) * (Ce*exp(-2*V/V_T) - C) / (exp(-2*V/V_T) - 1) : mole / second / meter**3
#electro_diffusion = -P_Ca * V / (R * V_T) * (Ce*exp(-2*V/V_T) - C) : mole / second / meter**3

C0 = Crest + (V0 - Vrest) * Cm / (F * s)   : mole / meter**3

# Figure out which branc to use: +sqrt or -sqrt? 
V0 = (B + sqrt(B**2 - 4.*A*CC)) / A        : volt   # at compartment edges

# Calcium diffusion
dC/dt = coupling + coupling_electro + electro_diffusion        : mole / meter**3
'''

N_astro = 2 # Total number of astrocytes in the network
astrocytes = NeuronGroup(N_astro, astro_eqs)

# Initital Conditions
astrocytes.L = 8 *  umeter  # length of a compartment
astrocytes.R = 3 * umeter
astrocytes.r = 2 * umeter
astrocytes.A = 1 / volt
astrocytes.C = 1.1e-4 * mole / meter**3

astro_mon = StateMonitor(astrocytes, variables=['C','coupling', 'coupling_electro', 'electro_diffusion', 'A', 'V0'], record=True)
b_mon = StateMonitor(astrocytes, variables=['VVT','V'], record=True)

# Diffusion between astrocytes
astro_to_astro_eqs = '''
A_pre = dR2_post / (L_post * s_post * V_T) : 1/volt (summed)
B_pre = 0.5 * (dR2_post / (L_post * s_post)) * (Vrest + V_post)/V_T : 1 (summed)
CC_pre = (dR2_post/L_post * s_post) * ((Crest-C_post) * (F * s_post/Cm) + Vrest * (V_post/V_T-1.)) : volt (summed)

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

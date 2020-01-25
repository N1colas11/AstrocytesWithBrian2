
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

duration = 500*ms          # Total simulation time
sim_dt = 10*ms            # Integrator/sampling step

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
#A                  : 1/volt
#B                  : 1
#CC                 : volt
A1                 : meter**4 / mole
B1                 : meter
C1                 : mole / meter**2
coupling_C         : mole / second / meter**3
coupling_electro   : mole / second / meter**3

V = Vrest + (F*s/Cm) * (C-Crest) : volt
VVT = -2.*V/V_T : 1  # value of 0, which leads to a zero denominator

electro_diffusion = -P_Ca * V / (R * V_T) * (Ce*exp(-2*V/V_T) - C) / (exp(-2*V/V_T) - 1) : mole / second / meter**3

# Figure out which branch to use: +sqrt or -sqrt? 
# or -sqrt... But C0 > 0 (is C1*A1 always greater than 0?)
C0 = (-B1 + sqrt(B1**2 + 4*C1*A1)) / (2.*A1) : mole / meter**3 
#C0 = B1 / A1 : mole / meter**3

V0 = Vrest + (C0 - Crest) * (F * s) / Cm   : volt

# Calcium diffusion
dC/dt = coupling_C + 0.*coupling_electro + 1.*electro_diffusion        : mole / meter**3
'''

N_astro = 2 # Total number of astrocytes in the network
astrocytes = NeuronGroup(N_astro, astro_eqs, method='euler')

# Initital Conditions
astrocytes.L = 8 *  umeter  # length of a compartment
astrocytes.R = 3 * umeter
astrocytes.r = 2 * umeter
astrocytes.C = [1.1e-4, 1.5e-4] * mole / meter**3

astro_mon = StateMonitor(astrocytes, variables=['C','coupling_C', 'coupling_electro', 'electro_diffusion', 'A1', 'B1', 'C1', 'V0'], record=True)
b_mon = StateMonitor(astrocytes, variables=['VVT','V', 'C0', 'dR2', 's', 'L'], record=True)

# Diffusion between astrocytes
astro_to_astro_eqs = '''
# Solution to quadratic voltage equation
#A_pre = dR2_post / (L_post * s_post * V_T) : 1/volt (summed)
#B_pre = 0.5 * (dR2_post / (L_post * s_post)) * (Vrest + V_post)/V_T : 1 (summed)
#CC_pre = (dR2_post/L_post * s_post) * ((Crest-C_post) * (F * s_post/Cm) + Vrest * (V_post/V_T-1.)) : volt (summed)

# Assuming the above is correct, let us figure out units. 
# using C_pre:  C * F * s / Cm = V  ==> F / Cm = V / s / C = (V / s) * m^3 / mole

# Solution to quadratic C0 equation (calcium concentration)
# A1 C0^2 + B1 C - C1 = 0
# C0 = (-B1 \pm \sqrt{B1^2 + 4*C1*A1) / (2*A1)
A1_pre = (0.5 * F * s_post)/(Cm * V_T) * (dR2_post / L_post) : meter**4 / mole (summed)
B1_pre = (1. - (s_post * F * C_post)/(2. * Cm * V_T)) : meter (summed)
C1_pre = dR2_post * C_post / L_post : mole / meter**2 (summed)

# Investigate why this should be a minus sign
coupling_C_post = (4*D_C / L_post**2) * (C_pre - C_post) : mole / second / meter**3 (summed)

coupling_electro_pre = (4*D_C/L_post**2/V_T) * (C0_post + C_post) * (V0_post - V_post) : mole/second/meter**3 (summed)
#coupling_pre = (4*D_C / L_post**2) * (C_pre - C_post) : mole / second / meter**3 (summed)
#coupling_electro_pre = (4*D_C/L_post**2/V_T) * (C0_post + C_post) * (V0_post - V_post) : mole/second/meter**3 (summed)
'''

astro_to_astro = Synapses(astrocytes, astrocytes, model=astro_to_astro_eqs, method='euler')
astro_to_astro.connect()

# Run Simulation
run(duration, report='text')

print("C= ", astro_mon.C)
print("C0= ", b_mon.C0)
print("coupling_C= ", astro_mon.coupling_C)
print("coupling_electro= ", astro_mon.coupling_electro)
print("electro_diffusion= ", astro_mon.electro_diffusion)
print("VVT=-2*V/V_T ", b_mon.VVT)
print("V0= ", astro_mon.V0)
print("V= ", b_mon.V)
print("V_T= ", V_T)
print("A1= ", astro_mon.A1)
print("B1= ", astro_mon.B1)
print("C1= ", astro_mon.C1)
print("s= ", b_mon.s)
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

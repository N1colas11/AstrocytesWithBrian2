# sent by Maurizio De Pitta on 2019-11-16

import numpy as np
from brian2 import *

import matplotlib.pylab as plt

# set_device('cpp_standalone', directory='./codegen', build_on_run=True)  # Use fast "C++ standalone mode"
# set_device('cython')
codegen.target = 'cython'

# A simple lambda just for parameter equivalence
rate_ = lambda x,L,r : x*L*(1-r)

################################################################################
# Model parameters
################################################################################
### General parameters
duration = 60*second       # Total simulation time
sim_dt = 50*ms               # Integrator/sampling step

### Astrocyte parameters
# ---  Calcium fluxes
# Value provided by MDP
O_P = 1.0*umolar/second      # Maximal Ca^2+ uptake rate by SERCAs  (0.9 in MDP, 1.0 in Evan)
# Value must equal v_er in Evan model
#O_P = 4.4*umolar/second      # Maximal Ca^2+ uptake rate by SERCAs  (0.9 in MDP, 1.0 in Evan)

K_P = 0.05 * umolar          # Ca2+ affinity of SERCAs
C_T = 2*umolar               # Total cell free Ca^2+ content
rho_A = 0.18                 # ER-to-cytoplasm volume ratio
Omega_C = 6/second           # Maximal rate of Ca^2+ release by IP_3Rs
Omega_L = 0.1/second         # Maximal rate of Ca^2+ leak from the ER
# --- IP_3R kinectics
d_1 = 0.13*umolar            # IP_3 binding affinity
d_2 = 1.05*umolar            # Ca^2+ inactivation dissociation constant
O_2 = 0.2/umolar/second      # IP_3R binding rate for Ca^2+ inhibition
d_3 = 0.9434*umolar          # IP_3 dissociation constant
d_5 = 0.08*umolar            # Ca^2+ activation dissociation constant
# --- IP_3 production
O_delta = 0.6*umolar/second  # Maximal rate of IP_3 production by PLCdelta
kappa_delta = 1.5* umolar    # Inhibition constant of PLC_delta by IP_3
K_delta = 0.1*umolar         # Ca^2+ affinity of PLCdelta
# --- IP_3 degradation
O_5P = 0.05/second       # Maximal rate of IP_3 degradation by IP-5P
K_5P = 10*umolar
K_D = 0.7*umolar             # Ca^2+ affinity of IP3-3K
K_3K = 1.0*umolar            # IP_3 affinity of IP_3-3K
O_3K = 4.5*umolar/second     # Maximal rate of IP_3 degradation by IP_3-3K
# --- IP_3 diffusion
#F = 0.09*umolar/second    # Maximal exogenous IP3 flow
I_Theta = 0.3*umolar         # Threshold gradient for IP_3 diffusion
omega_I = 0.05*umolar        # Scaling factor of diffusion

################################################################################
# Additional and modified parameters
################################################################################
Lambda = 2100*umeter**3

O_delta = 0.6 * umolar * Lambda * (1-rho_A) / second
O_3K    = 4.5 * umolar * Lambda * (1-rho_A) / second
O_5P    = 0.05 * umolar * Lambda * (1-rho_A) / second
F       = 0.1 / second
D       = 0.05 / second  # setting D to zero has no effect. Why?

################################################################################
# Model definition
################################################################################
defaultclock.dt = sim_dt     # Set the integration time

### Astrocytes
astro_eqs = '''
vol_coef = Lambda*(1-rho_A) : meter**3
dI/dt = (J_ex + J_delta - J_3K - J_5P + J_coupling)/vol_coef : mmolar
J_delta = 0*O_delta/(1 + I/kappa_delta) * C**2/(C**2 + K_delta**2) : mole/second
J_3K = 0*O_3K * C**4/(C**4 + K_D**4) * I/(I + K_3K)                : mole/second
J_5P = 0*O_5P*I/(I + K_5P)                                         : mole/second

# Exogenous stimulation (rectangular wave with period of 50s and duty factor 0.4)
stimulus = int((t % (50*second))<20*second)                      : 1
delta_I_bias = I - I_bias*stimulus                               : mmolar

# external input (called J_beta in Evan's thesis)
J_ex = 1*(-F*delta_I_bias*vol_coef)                                  : mole/second 

# Diffusion between astrocytes
J_coupling                                                       : mole/second

# Ca^2+-induced Ca^2+ release:
dC/dt = J_r + J_leak - J_pump                                : mmolar
dh/dt = (h_inf - h)/tau_h                                 : 1
J_r =  Omega_C * m_inf**3 * h**3 * (C_T - (1 + rho_A)*C) : mmolar/second
J_leak =  Omega_L * (C_T - (1 + rho_A)*C)                     : mmolar/second
J_pump =  O_P * C**2/(C**2 + K_P**2)                          : mmolar/second
m_inf = I/(I + d_1) * C/(C + d_5)                         : 1
h_inf = Q_2/(Q_2 + C)                                     : 1
tau_h = 1/(O_2 * (Q_2 + C))                               : second
Q_2 = d_2 * (I + d_1)/(I + d_3)                           : mmolar

# External IP_3 drive
I_bias : mmolar (constant)
'''

N_astro = 2 # Total number of astrocytes in the network
astrocytes = NeuronGroup(N_astro, astro_eqs, method='rk4')
astrocytes.I_bias = np.asarray([10, 0.],dtype=float)*umolar
astro_mon = StateMonitor(astrocytes, variables=['C', 'I', 'J_ex'], record=True)

# Diffusion between astrocytes
astro_to_astro_eqs = '''
coef = 0 : 1
delta_I = I_post - I_pre        : mmolar
J_coupling_post = -coef*D*delta_I*(Lambda*(1-rho_A))  : mole/second (summed)
'''
astro_to_astro = Synapses(astrocytes, astrocytes,
                          model=astro_to_astro_eqs)
syn_mon = StateMonitor(astro_to_astro, variables=['coef'], record=True)
astro_to_astro.connect(True)

################################################################################
# Simulation run
################################################################################
run(duration, report='text')

print(syn_mon.coef)
print("astro_to_astro= ", syn_mon.coef)

################################################################################
# Analysis and plotting
################################################################################
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(6.26894, 6.26894 * 0.66),
                       gridspec_kw={'left': 0.1, 'bottom': 0.12})
#ax.plot(astro_mon.t/second,astro_mon.C[0]/umolar,'b-')
#ax.plot(astro_mon.t/second,astro_mon.C[1]/umolar,'r-')
ax = axes[0]
ax.plot(astro_mon.t/second,astro_mon.C[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.C[1]/umolar,'r-')
ax.set_ylabel("CaCYT")
ax = axes[2]
ax.plot(astro_mon.t/second,astro_mon.I[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.I[1]/umolar,'r-')
ax.set_ylabel("I")
ax = axes[3]
ax.plot(astro_mon.t/second,astro_mon.J_ex[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.J_ex[1]/umolar,'r-')
ax.set_ylabel("J_ex")
plt.tight_layout()

#plt.show()
plt.savefig("plot_MDP.pdf")

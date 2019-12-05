#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:22:06 2019

@author: tnt19a
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 17:15:28 2019

@author: tnt19a
"""

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
duration = 5*10**5 *ms          # Total simulation time
sim_dt = 25*ms             # Integrator/sampling step
sim_dt = 12*ms             # Integrator/sampling step
sim_dt = 200*ms             # Integrator/sampling step

### Astrocyte parameters
# ---  Calcium fluxes
# Value provided by MDP
O_P = 1.0*umolar/second      # Maximal Ca^2+ uptake rate by SERCAs  (0.9 in MDP, 1.0 in Evan)
# Value must equal v_er in Evan model
#O_P = 4.4*umolar/second      # Maximal Ca^2+ uptake rate by SERCAs  (0.9 in MDP, 1.0 in Evan)

K_P = 0.05 * umolar          # Ca2+ affinity of SERCAs
C_T = 2*umolar               # Total cell free Ca^2+ content
rho_A = 0.18                 # ER-to-cytoplasm volume ratio
rho = rho_A / (1.+rho_A)     # ER-to-cytoplasm volume ratio
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
K_delta = 0.1*umolar         # Ca^2+ affinity of PLCdelta   <<<<<<<<<<<<<<<<<
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
# Volume of an average soma
Lambda = 2100*umeter**3

# Multiply value of amplitudes by cytosolic volume (changed rho_A to rho)
O_delta = 0.6  * umolar * Lambda * (1-rho) / second
O_3K    = 4.5  * umolar * Lambda * (1-rho) / second
O_5P    = 0.05 * umolar * Lambda * (1-rho) / second
F       = 0.1 / second
D       = 0.05 / second  # setting D to zero has no effect. Why?

print("**** Omega_L= ", Omega_L)
print("**** O_delta= ", O_delta) # 0.6
print("**** K_delta= ", K_delta) # 0.6

################################################################################
# Model definition
################################################################################
defaultclock.dt = sim_dt     # Set the integration time

### Astrocytes
astro_eqs = '''
vol_coef = Lambda*(1-rho) : meter**3
dI/dt = (J_ex + J_delta - J_3K - J_5P + J_coupling)/vol_coef : mmolar
J_delta = 1*O_delta/(1 + I/kappa_delta) * C**2/(C**2 + K_delta**2) : mole/second

J_3K = 1*O_3K * C**4/(C**4 + K_D**4) * I/(I + K_3K)                : mole/second  # ERROR? 
#J_3K = 0*O_3K * (k_d**4 / (CaCYT**4+k_d**4)) * (I / (I + k_3K))   : mole/second  # I not the same as MDP

J_5P = 1*O_5P*I/(I + K_5P)                                         : mole/second

# Exogenous stimulation (rectangular wave with period of 50s and duty factor 0.4)
# Pass stimulus through a smoothing function  exp(H(x-x0) - H(x-x1))
# period = 50*second
# length = 20*second
# 0->1 at t=0, 50, 100, ...
# 1->0 at t=20, 70, 120, ...
stimulus = int((t % (50*second))<20*second)                  : 1
#stimulus = sin(2*pi*t/(100*second))**10                               : 1
delta_I_bias = I - I_bias*stimulus                           : mmolar

# external input (called J_beta in Evan's thesis)
J_ex = 1*(-F*delta_I_bias*vol_coef)                          : mole/second 

# Diffusion between astrocytes
J_coupling                                                   : mole/second

# Ca^2+-induced Ca^2+ release:  (WHY IS CALCIUM ALWAYS ZERO? Currents Are not zero!!)
dC/dt = (J_chan + J_leak - J_pump)                           : mmolar
dh/dt = 1.*(h_inf - h)/tau_h                                 : 1

J_chan =  Omega_C * m_inf**3 * h**3 * (C_T - (1 + rho_A)*C)   : mmolar/second
J_leak =  1*Omega_L * (C_T - (1 + rho_A)*C)                     : mmolar/second
J_pump =  O_P * C**2/(C**2 + K_P**2)                          : mmolar/second

m_inf = I/(I + d_1) * C/(C + d_5)                            : 1
h_inf = Q_2/(Q_2 + C)                                        : 1
tau_h = 1/(O_2 * (Q_2 + C))                                  : second
Q_2 = d_2 * (I + d_1)/(I + d_3)                              : mmolar

# External IP_3 drive
I_bias : mmolar (constant)
'''

N_astro = 3 # Total number of astrocytes in the network
astrocytes = NeuronGroup(N_astro, astro_eqs, method='rk4')
astrocytes.I_bias = np.asarray([1, 0., 0.],dtype=float)*umolar
astro_mon = StateMonitor(astrocytes, variables=['C', 'J_chan', 'J_leak', 'J_pump', 'h',
                                                'I', 'J_delta', 'J_ex', 'J_3K', 'J_5P', 'J_coupling'], record=True)

# Diffusion between astrocytes
astro_to_astro_eqs = '''
coef = 1 : 1
delta_I = I_post - I_pre        : mmolar
J_coupling_post = -coef*D*delta_I*(Lambda*(1-rho_A))  : mole/second (summed)
'''
astro_to_astro = Synapses(astrocytes, astrocytes,
                          model=astro_to_astro_eqs)
syn_mon = StateMonitor(astro_to_astro, variables=['coef'], record=True)
# Branching Framework defined through connect values
#astro_to_astro.connect(True)

# connect compartmant soma=0<-->1 , 1<-->2
astro_to_astro.connect(i=0,j=1)
astro_to_astro.connect(i=1,j=2)

################################################################################
# Simulation run
################################################################################
run(duration, report='text')

#print("coef= ", syn_mon.coef)
print("C= ", astro_mon.C[0])

################################################################################
# Analysis and plotting
################################################################################
fig, axes = plt.subplots(nrows=N_astro, ncols=1, figsize=(6.26894, 8 * 6.26894 * 0.66),
                       gridspec_kw={'left': 0.1, 'bottom': 0.12})
for i in range(N_astro):
    ax = axes[i]
    ax.plot(astro_mon.t/second,astro_mon.C[i]/umolar)
    ax.set_ylabel("CaCYT")


plt.tight_layout()

plt.savefig("plot_s-1-1.pdf")

fig_in, axes_in = plt.subplots(nrows=N_astro, ncols=1, figsize=(6.26894, 8 * 6.26894 * 0.66),
                       gridspec_kw={'left': 0.1, 'bottom': 0.12})
for i in range(N_astro):
    ax = axes_in[i]
    ax.plot(astro_mon.t/second,astro_mon.C[i]/umolar)
    ax.set_ylabel("J_ex")


plt.tight_layout()

plt.savefig("plot_s-1-1_input.pdf")
plt.show()



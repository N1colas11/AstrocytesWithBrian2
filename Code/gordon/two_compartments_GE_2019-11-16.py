# sent by Maurizio De Pitta on 2019-11-16

# TODO
# Reinstate Diffusion
# Duplicate Evan's results from dissertation (comes first)
# Look into volume terms according to De Pitta's notes

import numpy as np
from brian2 import *

import matplotlib.pylab as plt

# set_device('cpp_standalone', directory='./codegen', build_on_run=True)  # Use fast "C++ standalone mode"
# set_device('cython')
codegen.target = 'cython'

# A simple lambda just for parameter equivalence
rate_ = lambda x,L,r : x*L*(1-r)

################################################################################
# Model parameters (consistent with C++ and dissertation)
################################################################################
### General parameters
duration = 30*second       # Total simulation time
sim_dt = 50*ms               # Integrator/sampling step

### Astrocyte parameters
# ---  Calcium fluxes
# Using my table from oschman_model.tex
r_c = v_IPR = Omega_C = 6/second
r_l = v_l = 0.11/second
v_er = 4.4*umolar/second
k_er = 0.05*umolar

d_1 = 0.13*umolar            # IP_3 binding affinity
d_2 = 1.049*umolar            # Ca^2+ inactivation dissociation constant
d_3 = 0.9434*umolar          # IP_3 dissociation constant
d_5 = 0.08234*umolar            # Ca^2+ activation dissociation constant

a_2 = 0.2/second/umolar     # Wrong in Evan (2019) and Maurizio (2009) and perhaps Li and Rinzel, 
                            #Appendix B (1994)

k_ER = 0.05*umolar
k_d = 0.7*umolar
k_delta = 1*umolar
K_delta = 0.5*umolar
k_pi = 0.6*umolar
k_3 = k_3K = 1.*umolar
k_R = 1.3 # units unclear
k_p = k_5P = 10. # units unclear

O_delta = 0.025*umolar/second
O_3K = 2.*umolar/second
O_5P = 0.05/second

Li = 0.4
Rb = 0.4
#rho_Ai = 

Dc = 0. * Hz
Dce = 0. * Hz
DI = 0. * Hz
D       = 0.05 / second

rho_A = 0.18  # revisit

"""
O_P = 0.9*umolar/second      # Maximal Ca^2+ uptake rate by SERCAs
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
F = 0.09*umolar/second    # Maximal exogenous IP3 flow
I_Theta = 0.3*umolar         # Threshold gradient for IP_3 diffusion
omega_I = 0.05*umolar        # Scaling factor of diffusion
"""

################################################################################
# Additional and modified parameters
################################################################################
Lambda = 2100*umeter**3

"""
O_delta = 0.6 * umolar * Lambda * (1-rho_A) / second
O_3K    = 4.5 * umolar * Lambda * (1-rho_A) / second
O_5P    = 0.05 * umolar * Lambda * (1-rho_A) / second
F       = 0.1 / second
D       = 0.05 / second
"""

################################################################################
# Model definition
################################################################################
defaultclock.dt = sim_dt     # Set the integration time

@check_units(x=mmolar, e=mmolar, n=1, result=1)
def Hill(x, e, n):
	return x**n / (x**n + e**n)
Hill.stateless = True

### Astrocytes
astro_eqs_evan_thesis = '''
vol_coef = Lambda*(1-rho_A) : meter**3
#rhoA = rho / (1-rho) : 1
dCaCYT/dt = Jchan + Jleak - Jpump : mmolar
dCaER/dt  = -(Jchan + Jleak - Jpump) / rho_A : mmolar
# Input J_ex or J_beta missing
dI/dt = (J_delta - J_3K - J_5P ) : mmolar   # return to J_coupling once code works
#dI/dt = (J_delta - J_3K - J_5P + J_coupling) : mmolar
dh/dt = (h_inf - h)/tau_h                                 : 1

J_delta = O_delta * Hill(k_delta, I, 1) * Hill(CaCYT, K_delta, 2.): mmolar/second
J_3K = O_3K * Hill(CaCYT, k_d, 4) * Hill(I, k_3K, 1.)             : mmolar/second  
J_5P = O_5P * I                                                   : mmolar/second

Q_2 = d_2 * (I + d_1)/(I + d_3)                           : mmolar
minf = Hill(I, d_1, 1)   : 1
ninf = Hill(CaCYT, d_5, 1)  : 1
Jchan = r_c * minf * ninf * h**3 * (CaER - CaCYT)     : mmolar/second
Jleak = r_l * (CaER - CaCYT)                          : mmolar/second
Jpump = v_er * Hill(CaCYT, k_er, 2)                   : mmolar/second

h_inf = Hill(Q_2, CaCYT, 1)                             : 1
tau_h = 1/(a_2 * (Q_2 + CaCYT))                         : second

# Exogenous stimulation (rectangular wave with period of 50s and duty factor 0.4)
stimulus = int((t % (50*second))<20*second)                       : 1
delta_I_bias = I - I_bias*stimulus                                : mmolar
# external input (called J_beta in Evan's thesis)
#J_ex = -F*delta_I_bias*vol_coef                                   : mole/second
# Diffusion between astrocytes
J_coupling                                                        : mole/second


# External IP_3 drive
I_bias : mmolar (constant)
'''

"""
astro_eqs_MDP = '''
vol_coef = Lambda*(1-rho_A) : meter**3
dI/dt = (J_delta - J_3K - J_5P + J_ex + J_coupling)/vol_coef : mmolar
J_delta = O_delta/(1 + I/kappa_delta) * C**2/(C**2 + K_delta**2) : mole/second
J_3K = O_3K * C**4/(C**4 + K_D**4) * I/(I + K_3K)                : mole/second
J_5P = O_5P*I/(I + K_5P)                                         : mole/second
# Exogenous stimulation (rectangular wave with period of 50s and duty factor 0.4)
stimulus = int((t % (50*second))<20*second)                      : 1
delta_I_bias = I - I_bias*stimulus                               : mmolar
# external input (called J_beta in Evan's thesis)
J_ex = -F*delta_I_bias*vol_coef                                  : mole/second
# Diffusion between astrocytes
J_coupling                                                       : mole/second

# Ca^2+-induced Ca^2+ release:
dC/dt = J_r + J_l - J_p                                   : mmolar
dh/dt = (h_inf - h)/tau_h                                 : 1
J_r = (Omega_C * m_inf**3 * h**3) * (C_T - (1 + rho_A)*C) : mmolar/second
J_l = Omega_L * (C_T - (1 + rho_A)*C)                     : mmolar/second
J_p = O_P * C**2/(C**2 + K_P**2)                          : mmolar/second
m_inf = I/(I + d_1) * C/(C + d_5)                         : 1
h_inf = Q_2/(Q_2 + C)                                     : 1
tau_h = 1/(O_2 * (Q_2 + C))                               : second
Q_2 = d_2 * (I + d_1)/(I + d_3)                           : mmolar

# External IP_3 drive
I_bias : mmolar (constant)
'''
"""

N_astro = 2 # Total number of astrocytes in the network
astrocytes = NeuronGroup(N_astro, astro_eqs_evan_thesis, method='rk4')
astrocytes.I_bias = np.asarray([10, 0.],dtype=float)*umolar
astro_mon = StateMonitor(astrocytes, variables=['CaCYT'], record=True)

astrocytes.CaCYT = [0.9,0.0]*umolar
astrocytes.CaER = [7.5,7.5]*umolar
astrocytes.h =  [0.95,0.95]
astrocytes.I = [0.100,0.2]*umolar  # umolar   # Units: umolar or ymolar? 

# Diffusion between astrocytes
astro_to_astro_eqs = '''
coef = 1 : 1
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
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6.26894, 6.26894 * 0.66),
                       gridspec_kw={'left': 0.1, 'bottom': 0.12})
ax.plot(astro_mon.t/second,astro_mon.CaCYT[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.CaCYT[1]/umolar,'r-')

plt.savefig("plot.pdf")

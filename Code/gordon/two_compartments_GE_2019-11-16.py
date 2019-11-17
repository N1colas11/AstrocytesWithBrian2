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
duration = 100*second       # Total simulation time
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

                            #Appendix B (1994)
k_ER = 0.05*umolar
k_d = 0.7*umolar
kappa_delta = 1*umolar
K_delta = 0.5*umolar
k_pi = 0.6*umolar
k_3 = k_3K = 1.*umolar
k_R = 1.3 # units unclear
k_p = k_5P = 10. * umolar  # units unclear

# Maximum amplitudes of current
#O_beta = 1.*umolar/second
O_delta = 0.025*umolar/second
O_3K = 2.*umolar/second
O_5P = 0.05/second
O_P = 1.*umolar/second   # J_SERCA (pump)
Omega_L = 1./second   # Leak frequency
Omega_C = 1./second   # IPR current frequency
Omega_2 = a_2 = 0.2/second/umolar     # Wrong units in Evan (2019) and Maurizio (2009) and perhaps Li and Rinzel, 
O_2 = 1.*umolar/second # 

# Morphological parameters
Li = 0.4
Rb = 0.4
#rho_Ai = 

# Diffusion coefficients
Dc = 0. * Hz
Dce = 0. * Hz
DI = 0. * Hz
D       = 0.05  * 0.001 / second
rho_A   = 0.18  # revisit

# IP3 external production
F = 0.2 * Hz                 # Maximal exogenous IP3 flow
I_Theta = 0.3*umolar         # Threshold gradient for IP_3 diffusion

################################################################################
# Additional and Modified Parameters
################################################################################
Lambda = 2100*umeter**3
# Different values than Evan in his thesis
O_delta = 0.6  * umolar * Lambda * (1-rho_A) / second  # DERIVED HOW?
O_3K    = 4.5  * umolar * Lambda * (1-rho_A) / second
O_5P    = 0.05 * umolar * Lambda * (1-rho_A) / second
F       = 0.1  * Hz   # Maximal exogeneous IP3 flow
D       = 0.05 / second

################################################################################
# Model definition
################################################################################
defaultclock.dt = sim_dt     # Set the integration time

# To go much faster, expand all Hill functions
# Alternatively, implement the Hill function within the C++ code generated. 
@check_units(x=mmolar, e=mmolar, n=1, result=1)
def Hill(x, e, n):
	return x**n / (x**n + e**n)
Hill.stateless = True

### Astrocytes
astro_eqs_evan_thesis = '''
# WHAT IS RHO? CONSTANT THROUGHOUT?
rho : 1
vol_coef = Lambda*(1-rho) : meter**3
rhoA = rho / (1-rho) : 1

dCaCYT/dt = Jchan + Jleak - Jpump : mmolar
dCaER/dt  = -(Jchan + Jleak - Jpump) / rhoA : mmolar
# Input J_ex or J_beta missing
dI/dt = (J_ex + J_delta - J_3K - J_5P ) / vol_coef : mmolar   # return to J_coupling once code works
#dI/dt = (J_delta - J_3K - J_5P + J_coupling) : mmolar
dh/dt = (h_inf - h) * Omega_h                                  : 1

#J_delta = O_delta * Hill(kappa_delta, I, 1) * Hill(CaCYT, K_delta, 2.): mole/second
J_delta = O_delta * (kappa_delta / (kappa_delta+I)) * (CaCYT**2 / (CaCYT**2+K_delta**2)) : mole/second

#J_3K = O_3K * Hill(k_d, CaCYT, 4) * Hill(I, k_3K, 1.)          : mole/second  
J_3K = O_3K * (k_d**4 / (CaCYT**4+k_d**4)) * (I / (I + k_3K))   : mole/second  

# Simplified version: J_5P = rbar * I (units of rbar are then Hz
#J_5P = O_5P * Hill(I, k_5P, 1)                                    : mole/second
J_5P = O_5P * (I / (I + k_5P))      : mole/second

#J_delta = O_delta/(1 + I/kappa_delta) * C**2/(C**2 + K_delta**2) : mole/second # MDP
#J_3K = O_3K * C**4/(C**4 + K_D**4) * I/(I + K_3K)                : mole/second # MDP
#J_5P = O_5P*I/(I + K_5P)                                         : mole/second # MDP

Q_2 = d_2 * (I + d_1)/(I + d_3)                           : mmolar
#minf = Hill(I, d_1, 1)   : 1
minf = I / (I + d_1) : 1

#ninf = Hill(CaCYT, d_5, 1)  : 1
ninf = CaCYT / (CaCYT + d_5) : 1
Jchan = r_c * (minf*ninf*h)**3 * (CaER - CaCYT)     : mmolar/second
Jleak = r_l * (CaER - CaCYT)                          : mmolar/second

#Jpump = v_er * Hill(CaCYT, k_er, 2)                   : mmolar/second
Jpump = v_er * CaCYT**2 / (CaCYT**2 + k_er**2)      : mmolar/second

#h_inf = Hill(Q_2, CaCYT, 1)                           : 1
h_inf = Q_2 / (Q_2 + CaCYT)                : 1

#Omega_h = Omega_2 * (Q_2 + CaCYT)                     : Hz
Omega_h = a_2 * (Q_2 + CaCYT)                        : Hz
#tau_h = 1/(a_2 * (Q_2 + CaCYT))                      : second  #(REVISIT CHOICES OF COEFFICIENTS)

# Exogenous stimulation (rectangular wave with period of 50s and duty factor 0.4)
stimulus = int((t % (50*second))<20*second)                       : 1
delta_I_bias = I - I_bias*stimulus                                : mmolar
# external input (called J_beta in Evan's thesis)
J_ex = -F * delta_I_bias * vol_coef                               : mole/second

# Diffusion between astrocytes
J_coupling                                                        : mole/second


# External IP_3 drive
I_bias : mmolar (constant)
'''

N_astro = 2 # Total number of astrocytes in the network
astrocytes = NeuronGroup(N_astro, astro_eqs_evan_thesis, method='rk4')
astrocytes.I_bias = np.asarray([10, 0.],dtype=float)*umolar
astro_mon = StateMonitor(astrocytes, variables=['CaCYT', 'I'], record=True)

astrocytes.CaCYT = [0.9,0.0]*umolar
#astrocytes.CaER = [7.5,7.5]*umolar
#astrocytes.h =  [0.95,0.95]
astrocytes.I = [0.100,0.2]*umolar  # umolar   # Units: umolar or ymolar? 
astrocytes.rho = 0.2
astrocytes.I_bias = np.asarray([10, 0.],dtype=float)*umolar

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
print("CaCYT: ", astro_mon.CaCYT)
print("I: ", astro_mon.I)

################################################################################
# Analysis and plotting
################################################################################
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6.26894, 6.26894 * 0.66),
                       gridspec_kw={'left': 0.1, 'bottom': 0.12})
ax = axes[0]
ax.plot(astro_mon.t/second,astro_mon.CaCYT[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.CaCYT[1]/umolar,'r-')
ax.set_ylabel("CaCYT")
ax = axes[1]
ax.plot(astro_mon.t/second,astro_mon.I[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.I[1]/umolar,'r-')
ax.set_ylabel("I")
plt.tight_layout()

plt.savefig("plot.pdf")

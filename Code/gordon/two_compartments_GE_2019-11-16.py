# sent by Maurizio De Pitta on 2019-11-16

# TODO
# Reinstate Diffusion
# Duplicate Evan's results from dissertation (comes first)
# Look into volume terms according to De Pitta's notes

import numpy as np
from brian2 import *

import matplotlib.pylab as plt

print("\n========== NEW SIMULATION: GE =========================================")

# set_device('cpp_standalone', directory='./codegen', build_on_run=True)  # Use fast "C++ standalone mode"
# set_device('cython')
codegen.target = 'cython'

# A simple lambda just for parameter equivalence
rate_ = lambda x,L,r : x*L*(1-r)

################################################################################
# Model parameters (consistent with C++ and dissertation)
################################################################################
### General parameters
duration = 60*second       # Total simulation time
duration = 150*ms    # Total simulation time
sim_dt = 50*ms               # Integrator/sampling step

### Astrocyte parameters
# ---  Calcium fluxes
# Using my table from oschman_model.tex

# Values in Evan Dissertation
v_er = 4.4*umolar/second   # SERCA pump max amplitude
r_c = v_IPR = Omega_C = 6/second
r_l = v_l = 0.11/second

# Parameters for MDP implementation (no equation for CE, closed system)
Omega_C = 6/second
Omega_L = 0.1/second

# Ratio of ER to cytosolic volume (assumed constant across cells)
rho_A   = 0.18  # revisit  Vol(ER) / Vol(cytosol)
rho = rho_A / (1 + rho_A)

# Values for r_c and r_l to use to get results consistent with MDP
# MDP assumes a closed system. has equations for CaER
r_c = Omega_C * rho_A
r_l = Omega_L * rho_A

# Value provided by MDP, taking into account that his equations are expressed in terms of C_T
# there is no equation for C_ER (no need if diffusion of CaCYT and C_ER are zero)
#r_c = v_IPR = Omega_C = 6/second
#r_l = v_l = 0.11/second
v_er = 1.0*umolar/second   # SERCA pump max amplitude

k_er = 0.05*umolar

d_1 = 0.13*umolar            # IP_3 binding affinity
d_2 = 1.049*umolar            # Ca^2+ inactivation dissociation constant
d_3 = 0.9434*umolar          # IP_3 dissociation constant
d_5 = 0.08234*umolar            # Ca^2+ activation dissociation constant

                            #Appendix B (1994)
k_ER = 0.05*umolar
k_d = 0.7*umolar
kappa_delta = 1*umolar
K_delta = 0.5*umolar    #<<<<<<<<<<<<<<<<<< (value used in Evan's thesis)
K_delta = 100*umolar    # (value used in MDP code) WHICH TO CHOOSE? 
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
#Omega_L = 1./second   # Leak frequency
#Omega_C = 1./second   # IPR current frequency
print("***** 1-rho= ", 1-rho)
Omega_2 = a_2 = 0.2/second/umolar     # Wrong units in Evan (2019) and Maurizio (2009) and perhaps Li and Rinzel, 
O_2 = 1.*umolar/second # 

# Morphological parameters
Li = 0.4
Rb = 0.4

# Diffusion coefficients
Dc = 0. * Hz
Dce = 0. * Hz
D_I = 0. * Hz
D       = 0.05 / second

# IP3 external production
F = 0.2 * Hz                 # Maximal exogenous IP3 flow
I_Theta = 0.3*umolar         # Threshold gradient for IP_3 diffusion

################################################################################
# Additional and Modified Parameters
################################################################################
Lambda = 2100*umeter**3
# Different values than Evan in his thesis!!! Where does 0.6, 4.5, 0.0 come from? 
# Values in Evan dissertation: O_delta= 0.025, O3k =2, r5p_bar=0.05
O_delta = 0.6  * umolar * Lambda * (1-rho) / second  # DERIVED HOW?
O_3K    = 4.5  * umolar * Lambda * (1-rho) / second
O_5P    = 0.05 * umolar * Lambda * (1-rho) / second
F       = 0.1  * Hz   # Maximal exogeneous IP3 flow
D       = 0.05 / second

print("******* r_l= ", r_l)
print("******* rho_A= ", rho_A)
print("******* Omega_L= ", Omega_L)
print("******* O_delta= ", O_delta)
print("******* K_delta= ", K_delta)

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
#################################
#   Morphology
################################
# Lambda is the reference volume of the compartment, taken as the soma volume
rho_A                     : 1  # Ratio of ER to CyT vol
rho = rho_A / (1+rho_A)   : 1  # Ratio of ER to Compartment vol
vol_coef = Lambda*(1-rho) : meter**3  # Cytosolic volume

##############################
#  Governing equations
##############################
dCaCYT/dt = Jchan + Jleak - Jpump            : mmolar
dCaER/dt  = -(Jchan + Jleak - Jpump) / rho_A : mmolar
dh/dt = 1*(h_inf - h) * Omega_h                : 1

# Input J_ex or J_beta missing
dI/dt = 1*(J_ex + J_delta - J_3K - J_5P + J_coupling) / vol_coef : mmolar

###############################
#   Currents / Fluxes
###############################

J_delta = 1*O_delta * (kappa_delta / (kappa_delta+I)) * 
                      (CaCYT**2 / (CaCYT**2+K_delta**2))             : mole/second # I close to the same
J_3K = 1*O_3K * (CaCYT**4 / (CaCYT**4+k_d**4)) * (I / (I + k_3K))    : mole/second  # I not the same as MDP

# Simplified version: J_5P = rbar * I 
# (units of rbar are then Hz)
J_5P = 1*O_5P * I / (I + k_5P)                                       : mole/second  # I exactly the same  as MDP

minf = I / (I + d_1) : 1
ninf = CaCYT / (CaCYT + d_5) : 1

Jchan = 1*r_c * (minf*ninf*h)**3 * (CaER - CaCYT)     : mmolar/second
#  J_leak =  1*Omega_L * (C_T - (1 + rho_A)*C)   ## MDP code
Jleak = 1*r_l * (CaER - CaCYT)                          : mmolar/second
Jpump = 1*v_er * CaCYT**2 / (CaCYT**2 + k_er**2)      : mmolar/second

Q_2 = d_2 * (I + d_1)/(I + d_3)                           : mmolar
h_inf = Q_2 / (Q_2 + CaCYT)                           : 1

Omega_h = a_2 * (Q_2 + CaCYT)                         : Hz
#tau_h = 1/(a_2 * (Q_2 + CaCYT))                      : second  #(REVISIT CHOICES OF COEFFICIENTS)

# Exogenous stimulation (rectangular wave with period of 50s and duty factor 0.4)
stimulus = int((t % (50*second))<20*second)                       : 1
delta_I_bias = I - I_bias*stimulus                                : mmolar
# external input (called J_beta in Evan's thesis)
J_ex = 1*(-F * delta_I_bias * vol_coef)                           : mole/second

# Diffusion between astrocytes
J_coupling                                                        : mole/second


# External IP_3 drive
I_bias : mmolar (constant)
'''

N_astro = 2 # Total number of astrocytes in the network
astrocytes = NeuronGroup(N_astro, astro_eqs_evan_thesis, method='rk4')
astro_mon = StateMonitor(astrocytes, variables=['CaER', 'CaCYT', 'I', 'J_coupling', 'vol_coef',
   'J_ex', 'Jchan', 'Jleak', 'Jpump', 'rho_A', 'rho', 'h', 'J_delta', 'J_3K', 'J_5P'], record=True)

#######################################
#   Initial conditions
#######################################
#astrocytes.CaCYT = [0.9,0.0]*umolar
#astrocytes.CaER = [7.5,7.5]*umolar
#astrocytes.h =  [0.95,0.95]
#astrocytes.I = [0.100,0.2]*umolar  # umolar   # Units: umolar or ymolar? 
astrocytes.rho_A  = 0.18
astrocytes.I_bias = np.asarray([10, 0.],dtype=float)*umolar
astrocytes.CaER = 2.*umolar/rho_A   # C_T in MDP code

##########################################3
# Diffusion between astrocytes
##########################################3
astro_to_astro_eqs = '''
coef = 0 : 1
delta_I = I_post - I_pre        : mmolar
J_coupling_post = -coef * D * delta_I * vol_coef  : mole/second (summed)
'''
astro_to_astro = Synapses(astrocytes, astrocytes,
                          model=astro_to_astro_eqs)
syn_mon = StateMonitor(astro_to_astro, variables=['coef'], record=True)
astro_to_astro.connect(True)

################################################################################
# Simulation run
################################################################################
run(duration, report='text')

print("CaCYT: ", astro_mon.CaCYT[0])
print("CaER: ", astro_mon.CaER[0])
print("I: ", astro_mon.I[0])
print("Jchan= ", astro_mon.Jchan[0])
print("Jleak= ", astro_mon.Jleak[0])
print("Jpump= ", astro_mon.Jpump[0])
print("J_coupling: ", astro_mon.J_coupling[0])
print("vol_coef= ", astro_mon.vol_coef[0])
print("h= ", astro_mon.h[0])
print("J_delta= ", astro_mon.J_delta[0])
print("J_3K= ", astro_mon.J_3K[0])
print("J_5P= ", astro_mon.J_5P[0])
print("rho_A= ", rho_A)
print("monitor rho_A= ", astro_mon.rho_A[0])
print("rho= ", astro_mon.rho[0])

################################################################################
# Analysis and plotting
################################################################################
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(6.26894, 6.26894 * 0.66),
                       gridspec_kw={'left': 0.1, 'bottom': 0.12})
ax = axes[0]
ax.plot(astro_mon.t/second,astro_mon.CaCYT[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.CaCYT[1]/umolar,'r-')
ax.set_ylabel("CaCYT")
ax = axes[1]
ax.plot(astro_mon.t/second,astro_mon.CaER[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.CaER[1]/umolar,'r-')
ax.set_ylabel("CaER")
ax = axes[2]
ax.plot(astro_mon.t/second,astro_mon.I[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.I[1]/umolar,'r-')
ax.set_ylabel("I")
ax = axes[3]
ax.plot(astro_mon.t/second,astro_mon.J_ex[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.J_ex[1]/umolar,'r-')
ax.set_ylabel("J_ex")
#plt.tight_layout()

plt.savefig("plot.pdf", bbox_inches='tight')

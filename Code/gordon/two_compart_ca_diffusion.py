#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 18:28:40 2019

@author: tnt19a
"""

"""
Figure 5: Astrocytes connected in networks

De Pitta model with multiple compartments, but without diffusion. 
All compartments are the same, but there is a soma. 
"""
import matplotlib      # DELETE
# matplotlib.use('agg')  # DELETE
from brian2 import *
from brian2.units import *
from brian2.units.allunits import ymole, ymol, ymolar, umolar
#import brian2.units as un
#from brian2.units.constants import * # GE

# import plot_utils as pu

set_device('cpp_standalone', directory='./codegen', build_on_run=True)  # Use fast "C++ standalone mode"

################################################################################
# Model parameters
################################################################################
### General parameters
## Something wrong with the time scales. I should be able to have much larger time steps
duration = 4500*msecond       # Total simulation time
sim_dt = 40.*ms               # Integrator/sampling step (too small!!)  (beyond 0.30*ms, there is 2-point oscillation in solution)
#duration = 50*second       # Total simulation time
#sim_dt = 1*second               # Integrator/sampling step (too small!!)  (beyond 0.30*ms, there is 2-point oscillation in solution)
# one pico = (1000)^4 yocto

rho=0.5
Lambda=2100*umeter**3
Jbeta=0.1e-24*mole/second   # IP3 input
# kappad=1*umolar # <<<<< ???? --> kappa_delta

### Astrocyte parameters
# ---  Calcium fluxes
O_P = 4.4*umolar/second      # Maximal Ca^2+ uptake rate by SERCAs (Op) (0.9)
K_P = 0.05 * umolar          # Ca2+ affinity of SERCAs (Kp)
Omega_C = 6./second           # Maximal rate of Ca^2+ release by IP_3Rs
Omega_L = 0.11/second         # Maximal rate of Ca^2+ leak from the ER (0.1)

#-------------------------------
# Switches to turn diffusion on or off
cDCa = 1.
cDIP3 = 1.
# --- Calcium diffusion
D_Ca = cDCa * 1*Hz
# --- IP3 diffusion
D_IP3 = cDIP3 * 1*Hz     # TEMPORARY value
#-------------------------------

# --- IP_3R kinectics
d_1 = 0.13*mole #umolar            # IP_3 binding affinity
d_2 = 1.049*umolar            # Ca^2+ inactivation dissociation constant (1.05)
O_2 = 0.1335/umolar/second      # IP_3R binding rate for Ca^2+ inhibition (0.2)
d_3 = 0.9434*mole #umolar          # IP_3 dissociation constant
d_5 = 0.082*umolar            # Ca^2+ activation dissociation constant (0.08)
# --- IP_3 production
O_delta = 0.02e-24*mole/second  # Maximal rate of IP_3 production by PLCdelta  (0.06)
kappa_delta = 1.0*mole # umolar    # Inhibition constant of PLC_delta by IP_3  (kappad?) (consistent with ode file)
#kappa_delta = 1.5* umolar    # Inhibition constant of PLC_delta by IP_3  (kappad?)
K_delta = 0.1*umolar         # Ca^2+ affinity of PLCdelta
# --- IP_3 degradation
O_5P = 0.05e-24*mole/second       # Maximal rate of IP_3 degradation by IP-5P (O5p)
K_5P = 10.*mole #umolar
K_D = 0.7*umolar             # Ca^2+ affinity of IP3-3K
K_3K = 1.0*mole # umolar            # IP_3 affinity of IP_3-3K
O_3K = 2.0e-24*mole/second     # Maximal rate of IP_3 degradation by IP_3-3K  (4.5)
# --- IP_3 diffusion

#-------------------------------------------------------
# ODEs
#-------------------------------------------------------
# IP3   : cytosolic IP3
# CaCYT : cytosolic calcium
# CaER  : endoplasmic reticulum calcium
# h     : IP3R deinactivation
#-------------------------------------------------------

################################################################################
# Model definition
################################################################################
defaultclock.dt = sim_dt     # Set the integration time

# (used to multiply ip3 equations)  *1e3/(Lambda*(1-rho))  : mmolar
comp_eqs = '''
dIP3/dt=(Jbeta+Jdelta-J3k-J5p+Jdiff_IP3) : mole  # check sign of diffusion
dCaCYT/dt = Jcicr+Jleak-Jserca+Jdiff_Ca  : mmolar
dCaER/dt=(-Jcicr-Jleak+Jserca) / (rho / (1-rho) ) : mmolar
dh/dt=(hinf-h)/tauh : 1

#IP3 dynamics
Jdelta=O_delta * ( 1 - ( IP3 /(kappa_delta + IP3) ) * ( CaCYT**2/(CaCYT**2+K_delta**2) )  ) : mole/second
J3k=O_3K * ( ( CaCYT**4/(CaCYT**4+K_D**4)) * (IP3/(K_3K+IP3) ) ) : mole/second
J5p=O_5P * IP3 /(IP3 + K_5P) : mole/second

# Calcium core
minf=IP3/(IP3+d_1) : 1
ninf=CaCYT/(CaCYT+d_5) : 1
Jcicr=Omega_C*(minf*ninf*h)**3*(CaER-CaCYT) : mmolar/second
Jleak=Omega_L*(CaER-CaCYT) : mmolar/second
Jserca=O_P*CaCYT**2/(K_P**2+CaCYT**2) : mmolar/second

Jdiff_Ca : mmolar/second
Jdiff_IP3 : mole/second

# Gating variable
Q2=d_2*(IP3+d_1)/(IP3+d_3) : mmolar
hinf=Q2/(Q2+CaCYT) : 1
tauh=1./O_2/(Q2+CaCYT) : second
'''

#-------------------------------------------------------

N_comp = 2 # Total number of compartments in the network
compartments = NeuronGroup(N_comp, comp_eqs, method='rk4')

# Initial conditions
# How are these set more generally with random numgers? 
compartments.CaCYT = [0.9,0.0]*umolar
compartments.CaER = [7.5,7.5]*umolar
compartments.h =  [0.95,0.95]
compartments.IP3 = [0.100,0.0]*mole  # umolar   # Units: umolar or ymolar? 

# # Diffusion between compartments
synapse_eqs = '''
coef= 1 : 1
Jdiff_Ca_post = -coef * D_Ca * (CaCYT_post - CaCYT_pre) : mmolar/second (summed)
Jdiff_IP3_post = -D_IP3 * (IP3_post - IP3_pre) : mole/second (summed) 
# multiply by the volume of postsynatpic compartment: Lambda (1-rhoA)
'''


comp_to_comp = Synapses(compartments, compartments,
                        model=synapse_eqs)

comp_to_comp.connect()

################################################################################
# Monitors
################################################################################
astro_Ca_mon = StateMonitor(compartments, variables=['CaCYT'], record=True)
astro_IP3_mon = StateMonitor(compartments, variables=['IP3'], record=True)
astro_diff_mon = StateMonitor(compartments, variables=['Jdiff_Ca'], record=True)

syn_mon = StateMonitor(compartments, variables=['Jdiff_IP3'], record=True)


################################################################################
# Simulation run
################################################################################
run(duration, report='text')
#print(astro_diff_mon.Jdiff_Ca)
print("CaCYT, ", astro_Ca_mon.CaCYT)
print("IP3, ", astro_IP3_mon.IP3)
################################################################################
# Analysis and plotting
################################################################################
fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, figsize=(6.26894 * 2, 6.26894 * 0.66),
                       gridspec_kw={'left': 0.1, 'bottom': 0.12})

ax1.plot(astro_Ca_mon.t/second,astro_Ca_mon.CaCYT[0]/umolar,label='C1')
ax1.plot(astro_Ca_mon.t/second,astro_Ca_mon.CaCYT[1]/umolar,label='C2')
ax1.legend()
ax1.set(xlabel='time (ms)',ylabel='Ca concentration (umolar)')

ax2.plot(astro_IP3_mon.t/second,astro_IP3_mon.IP3[0]/umolar,label='C1')
ax2.plot(astro_IP3_mon.t/second,astro_IP3_mon.IP3[1]/umolar,label='C2')
ax2.legend()
ax2.set(xlabel='time (ms)',ylabel='IP3 concentration (umolar)')
plt.savefig("plot.pdf")
#plt.show()
